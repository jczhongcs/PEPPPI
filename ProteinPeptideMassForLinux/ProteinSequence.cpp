#include "stdafx.h"
#include "ProteinSequence.h"
#include <fstream>


#include <sstream>
#include <map>
#include <algorithm>
#include "PeptideTools.h"
#include <vector>

CProteinSequence::CProteinSequence(void)
{
}


CProteinSequence::~CProteinSequence(void)
{
}




ifstream& CProteinSequence::getcleanline( ifstream& is, string& buffer )
{
	if (getline(is,buffer))
	{
		if (buffer.size()>0&&buffer.at(buffer.size()-1)=='\r')
		{
			buffer.erase(buffer.size()-1);
		}
	}
	return is;
}

void CProteinSequence::parse_file( string strProteinSequenceName )
{
	ofstream outf13("proteinMass_Output.txt");

//	CMassPeptideIons* pMassPeptideIons=NULL;
	ifstream is (strProteinSequenceName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strProteinSequenceName).c_str());
	}
	int irow=0;
	int iPeptideRow=0;
	string iIsoformNumber;
	map< string, vector<string> > mapvPeptide;
	string strNameOld="";
	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{

		if (buffer[0]=='>')
		{
			istringstream strStrem(buffer);//×Ö·û´®±ä³ÉÁ÷
			string strName;
			strStrem>>strName;
			//strName=buffer;
			m_vstrName.push_back(strName);
			
			string strIsoName;
			strStrem>>strIsoName;


			string::size_type pos2 = strIsoName.find('_');
			iIsoformNumber=strIsoName.substr(pos2+1, strIsoName.length());
			if (strNameOld!=strName)
			{
				if (irow>0)//consider each protein has one or more isoform
				{
					m_vmapvIsoformPeptide.push_back(mapvPeptide);//add the isoforms of the protein 
				}
				strNameOld=strName;
				mapvPeptide.clear();
			}
			continue;
		}

		
		istringstream strStrem(buffer);//×Ö·û´®±ä³ÉÁ÷
		string strPeptide;
		vector<string> vstrPeptide;
		while (strStrem>>strPeptide)
		{
			if (strPeptide.length()>0)
			{
				vstrPeptide.push_back(strPeptide);
			}
		}
		mapvPeptide.insert(make_pair(iIsoformNumber,vstrPeptide));
		
		irow++;//protein no
		iPeptideRow++;
	}
	m_vmapvIsoformPeptide.push_back(mapvPeptide);//add last  isoforms of the protein
}

void CProteinSequence::parse_fileLine( string strProteinSequenceName )
{
	ofstream outf13("proteinMass_Output.txt");

	//	CMassPeptideIons* pMassPeptideIons=NULL;
	ifstream is (strProteinSequenceName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strProteinSequenceName).c_str());
	}

	string iIsoformNumber;
	string strNameOld="";
	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{

		if (buffer[0]=='>')
		{
			istringstream strStrem(buffer);//×Ö·û´®±ä³ÉÁ÷
			string strName;
			strStrem>>strName;
			//strName=buffer;
			m_vstrName.push_back(strName);

			string strIsoName;
			strStrem>>strIsoName;


			string::size_type pos2 = strIsoName.find('_');
			iIsoformNumber=strIsoName.substr(pos2+1, strIsoName.length());
			int b=atoi(iIsoformNumber.c_str());

			m_viIsoform.push_back(b);

			continue;
		}

		m_vstrProteinLine.push_back(buffer);
		
	}
}

int CProteinSequence::CountEnzyme(string &strPattern)
{
	int iCountEnzyme=0;
	for (int i=1;i<(strPattern.length()-1);i++)
	{
		if ((strPattern.at(i)=='R') || (strPattern.at(i)=='K'))
			iCountEnzyme++;
	}
	return iCountEnzyme;
}
int CProteinSequence::Alignment( string& strSource, string& strPattern)
{
	//judge the pattern has contain R K
	CPeptideTools objTools;
	if (strSource==strPattern)
	{
		return 1;
	}else if (objTools.SUNDAYFindSigle(const_cast<char*>(strSource.c_str()),const_cast<char*>(strPattern.c_str()))>0)
	{
		return 2;
	}else if (objTools.ldistance(strSource,strPattern)<2)
	{
		return 3;
	}else
	{
		return 0;
	}
	

	
/*

	objTools.massQuantity("LVTSGAESGNLNTSPSSNQTR");
	int i;
	i=objTools.ldistance("LVTSGAESGNLNTSPSSNQTR","LVTSGAESGNLNTSPSSNQTR");
	i=objTools.ldistance("LVTSGAESGNLNTSPSSNQTR","TSGAESGNLNTSPSSNQTR");
	i=objTools.ldistance("LVTSGAESGNLNTSPSSNQT","LVTSGAESGNLNTSPSSNQTR");
	i=objTools.ldistance("LVTSGAESGNLNSPSSNQTR","LVTSGAESGNLNTSPSSNQTR");
	i=objTools.KMP("LVTSGAESGNL","TSGAES");
	i=objTools.SUNDAYFindSigle("LVTSGAESGNL","TSGAESR");
*/

	return 0;

}

void CProteinSequence::PeptideProteinAlignment( CListMassPeptideXcorr* pListMassPeptideXcorr )
{
	ofstream outf13("proteinMassPeptide_Output.txt");


	for (int i=0;i<pListMassPeptideXcorr->m_vpMassPeptideXcorr.size();i++)//for each peptide
	{
		outf13<<i<<'\t'<<pListMassPeptideXcorr->m_vpMassPeptideXcorr[i]->m_strPeptideSequence<<endl;
		for (int j=0;j<m_vmapvIsoformPeptide.size();j++)//protein
		{
			
			//cout<<j<<endl;
			
			for (int k=1;k<=(m_vmapvIsoformPeptide[j].size());k++)//map, protein has many isforms
			{
				//int to stirng
				stringstream ss;
				ss<<k;
				string   s=ss.str();
				//////////////////////////////////////////////////////////////////////////
				//alignment peptide with protein peptide
				int iProteinPeptideSplitCount=m_vmapvIsoformPeptide[j][s].size();
				for (int l=0;l<iProteinPeptideSplitCount;l++)//isform has many peptide
				{
					
					string& strMassPeptide=pListMassPeptideXcorr->m_vpMassPeptideXcorr[i]->m_strPeptideSequence;
					string& strProtein=m_vmapvIsoformPeptide[j][s][l];
					string strProteinSequence(strProtein);
					int iCountEnzyme=CountEnzyme(strMassPeptide);
					for (int h=0;h<iCountEnzyme && l+h+1<iProteinPeptideSplitCount;h++) 
						strProteinSequence+=m_vmapvIsoformPeptide[j][s][l+h+1];//connect
					int iCount=strMassPeptide.length()-strProteinSequence.length();
					if (iCount>0) continue;
					int iAlignment=Alignment(strProteinSequence,strMassPeptide);
					if (iAlignment>0)
					{
						outf13<<i<<'\t'<<j<<'\t'<<k<<'\t'<<l<<'\t'<<iCountEnzyme<<'\t'<<strMassPeptide<<'\t'<<strProtein<<'\t'<<strProteinSequence<<'\t'<<iAlignment<<endl;
					}
					
				}
			}
			 

		}
	}
}


/*

void CProteinSequence::PeptideProteinAlignmentWithLine( CListMassPeptideXcorr* pListMassPeptideXcorr )
{
	ofstream outf13("proteinMassPeptideLine_Output.txt");
	CPeptideTools objTools;

	for (int i=0;i<pListMassPeptideXcorr->m_vpMassPeptideXcorr.size();i++)//for each peptide
	{

		int iName=0;

		outf13<<i<<'\t'<<pListMassPeptideXcorr->m_vpMassPeptideXcorr[i]->m_strPeptideSequence<<endl;
		for (int j=0;j<m_vmapvIsoformPeptide.size();j++)//protein
		{

			//cout<<j<<endl;

			
			for (map < string, vector<string> >::iterator m1_Iter = m_vmapvIsoformPeptide[j].begin( ); m1_Iter != m_vmapvIsoformPeptide[j].end( ); m1_Iter++ )
//			for (int k=1;k<=m_vmapvIsoformPeptide[j].size();k++)//map, protein has many isforms
			{
				//int to stirng
// 				stringstream ss;
// 				ss<<k;
// 				string   s=ss.str();
				//////////////////////////////////////////////////////////////////////////
				//alignment peptide with protein peptide

				
// 				int iProteinPeptideSplitCount=m1_Iter->second.size();
// 				for (int l=0;l<iProteinPeptideSplitCount;l++)//isform has many peptide
// 				{
// 

					string strMassPeptide=pListMassPeptideXcorr->m_vpMassPeptideXcorr[i]->m_strPeptideSequence;
					string& strProtein=m1_Iter->second[0];
					string strProteinSequence(strProtein);
					int iCount=strMassPeptide.length()-strProteinSequence.length();
					if (iCount>1) continue;
					std::transform(strMassPeptide.begin(), strMassPeptide.end(),strMassPeptide.begin(), ::toupper);//upper
//					int iAlignment=Alignment(strProteinSequence,strMassPeptide);
					objTools.SUNDAYFind(const_cast<char*>(strProteinSequence.c_str()),const_cast<char*>(strMassPeptide.c_str()),j,m1_Iter->first,outf13,m_vstrName[iName]);
					iName++;
// 					if (iAlignment>0)
// 					{
// 						outf13<<i<<'\t'<<j<<'\t'<<k<<'\t'<<l<<'\t'<<iCountEnzyme<<'\t'<<strMassPeptide<<'\t'<<strProtein<<'\t'<<strProteinSequence<<'\t'<<iAlignment<<endl;
// 					}

//				}
				
			}


		}
	}
}*/



void CProteinSequence::PeptideProteinAlignmentWithLine( CListMassPeptideXcorr* pListMassPeptideXcorr )
{
	ofstream outf13("proteinMassPeptideLine_Output.txt");
	CPeptideTools objTools;

	for (int i=0;i<pListMassPeptideXcorr->m_vpMassPeptideXcorr.size();i++)//for each peptide
	{

		
		string strMassPeptide=pListMassPeptideXcorr->m_vpMassPeptideXcorr[i]->m_strPeptideSequence;


		outf13<<i<<'\t'<<pListMassPeptideXcorr->m_vpMassPeptideXcorr[i]->m_strPeptideSequence<<endl;
		for (int j=0;j<m_vstrProteinLine.size();j++)//protein
		{

			string strProteinSequence(m_vstrProteinLine[j]);
			int iCount=strMassPeptide.length()-strProteinSequence.length();
			if (iCount>1) continue;
			std::transform(strMassPeptide.begin(), strMassPeptide.end(),strMassPeptide.begin(), ::toupper);//upper

			objTools.SUNDAYFind(const_cast<char*>(strProteinSequence.c_str()),const_cast<char*>(strMassPeptide.c_str()),j,m_viIsoform[j],outf13,m_vstrName[j]);
				
			


		}
	}


}

void CProteinSequence::parse_filePeptideProteinMapp( string strfilePeptideProteinMapp )
{

	ofstream outf13("proteinMass_Output.txt");

	//	CMassPeptideIons* pMassPeptideIons=NULL;
	ifstream is (strfilePeptideProteinMapp.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strfilePeptideProteinMapp).c_str());
	}
	int irow=0;
	int iPeptideRow=0;
	string iIsoformNumber;

	string strNameOld="";
	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
		if (irow++==0) continue;

		istringstream strStrem(buffer);

		string strProtienNames;
		strStrem>>strProtienNames;

		string strSpectrumName;
		strStrem>>strSpectrumName;

		string dXcorScore;
		strStrem>>dXcorScore;

		string strPeptideSequence;
		strStrem>>strPeptideSequence;

		string dPeptideMass;
		strStrem>>dPeptideMass;

		string iCharge;
		strStrem>>iCharge;

		string::size_type pos2 = strProtienNames.find(',');
		string strProteinNameRest=strProtienNames;
		while (pos2!=-1)
		{
			string strSingleProtein=strProteinNameRest.substr(0, pos2);
			outf13<<strSingleProtein<<'\t'<<strSpectrumName<<'\t'<<dXcorScore<<'\t'<<strPeptideSequence<<'\t'<<dPeptideMass<<'\t'<<iCharge<<endl;


			strProteinNameRest=strProteinNameRest.substr(pos2+1, strProteinNameRest.length());
			pos2 = strProteinNameRest.find(',');

		}

		outf13<<strProteinNameRest<<'\t'<<strSpectrumName<<'\t'<<dXcorScore<<'\t'<<strPeptideSequence<<'\t'<<dPeptideMass<<'\t'<<iCharge<<endl;


	}

		

}

void CProteinSequence::parse_fileProteinIPI_Swissport_Map( string strfilePeptideProteinMapp )
{

	ofstream outf13("parse_fileProteinIPI_Swissport_Map_output2008.txt");

	//	CMassPeptideIons* pMassPeptideIons=NULL;
	ifstream is (strfilePeptideProteinMapp.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strfilePeptideProteinMapp).c_str());
	}
	int irow=0;
	int iPeptideRow=0;
	string iIsoformNumber;

	string strNameOld="";
	string strIPIID="";
	string strSwissID="";
	string strGeneNam="";
	string strRefseq="";
	vector<string> vstrID;
	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
//		if (irow++==0) continue;

	


		if (buffer[0]=='I'&&buffer[1]=='D')
		{
			strIPIID=buffer;
			continue;
		}
		

		string::size_type pos1 = buffer.find("DR   Entrez Gene");

		if (pos1!=-1)
		{
			strGeneNam=strGeneNam+';'+buffer;irow++;
			continue;
		}

		string::size_type pos12 = buffer.find("DR   REFSEQ_REVIEWED");

		if (pos12!=-1)
		{
			strRefseq=strRefseq+';'+buffer;
			continue;
		}
		
		
		string::size_type pos2 = buffer.find("UniProtKB/Swiss-Prot");
		
		if (pos2!=-1)
		{

//			strSwissID=strSwissID+'\t'+buffer;
			vstrID.push_back(buffer);
			continue;
		}

		pos2 = buffer.find("UniProtKB/TrEMBL");

		if (pos2!=-1)
		{
//			strSwissID=strSwissID+'\t'+buffer;
			vstrID.push_back(buffer);
			
			continue;
		}


		if (buffer[0]=='/'&&buffer[1]=='/')
		{
			for (vector<string>::iterator it=vstrID.begin();
				it!=vstrID.end();
				it++)
			{
				outf13<<strIPIID<<'\t'<<strGeneNam<<'\t'<<strRefseq<<'\t'<<*it<<endl;
				
			}
			if (vstrID.size()==0)
			{
				outf13<<strIPIID<<'\t'<<strGeneNam<<'\t'<<strRefseq<<'\t'<<endl;
			}


			strIPIID="";
//			strSwissID="";
			strGeneNam="";
			strRefseq="";
			vstrID.clear();
			
		}
		
	}
	cout<<irow;
	cin>>irow;
}

void CProteinSequence::parse_filePPIMappingProteinIPI_Swissport_Map( string strfilePPIMapp )
{

	ofstream outf13("parse_filePPIMappingProteinIPI_Swissport_Map.txt");

	//	CMassPeptideIons* pMassPeptideIons=NULL;
	ifstream is (strfilePPIMapp.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strfilePPIMapp).c_str());
	}
	int irow=0;
	int iPeptideRow=0;


	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
		

		istringstream strStrem(buffer);

		string strProtienPPIID;
		strStrem>>strProtienPPIID;

		string strProteinName;
		strStrem>>strProteinName;

		string strNMID;
		strStrem>>strNMID;

		string strNPID;
		strStrem>>strNPID;

		string strentrezgeneID;
		strStrem>>strentrezgeneID;

		string stromim_id;
		strStrem>>stromim_id;

		string strswissprot_id;
		strStrem>>strswissprot_id;

		string strmain_name;
		strStrem>>strmain_name;

		string::size_type pos2 = strswissprot_id.find(',');
		string strSwissProtIDRest=strswissprot_id;
		while (pos2!=-1)
		{
			string strSingleProtein=strSwissProtIDRest.substr(0, pos2);
			outf13<<strProtienPPIID<<'\t'<<strProteinName<<'\t'<<strNMID<<'\t'<<strNPID<<'\t'<<strentrezgeneID<<'\t'<<stromim_id<<'\t'<<strSingleProtein<<'\t'<<strmain_name<<endl;


			strSwissProtIDRest=strSwissProtIDRest.substr(pos2+1, strSwissProtIDRest.length());
			pos2 = strSwissProtIDRest.find(',');

		}

		
		outf13<<strProtienPPIID<<'\t'<<strProteinName<<'\t'<<strNMID<<'\t'<<strNPID<<'\t'<<strentrezgeneID<<'\t'<<stromim_id<<'\t'<<strSwissProtIDRest<<'\t'<<strmain_name<<endl;


	}
	cout<<irow;
	cin>>irow;
}


void CProteinSequence::parse_PeptideProteinSequest( string strfilePPIMapp )
{

	ofstream outf13("parse_PPSequestData02Gt05.txt");

	//	CMassPeptideIons* pMassPeptideIons=NULL;
	ifstream is (strfilePPIMapp.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strfilePPIMapp).c_str());
	}
	int irow=0;
	int iPeptideRow=0;


	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{


		istringstream strStrem(buffer);

		string strProtienIPIID;
		strStrem>>strProtienIPIID;

		string strNumber;
		strStrem>>strNumber;


		string::size_type pos2 = strProtienIPIID.find(',');
		string strSwissProtIDRest=strProtienIPIID;
		while (pos2!=-1)
		{
			string strSingleProtein=strSwissProtIDRest.substr(0, pos2);
			outf13<<strSingleProtein<<'\t'<<strNumber<<endl;


			strSwissProtIDRest=strSwissProtIDRest.substr(pos2+1, strSwissProtIDRest.length());
			pos2 = strSwissProtIDRest.find(',');

		}


		outf13<<strSwissProtIDRest<<'\t'<<strNumber<<endl;


	}
	cout<<irow;
//	cin>>irow;
}


void CProteinSequence::parse_ProteinIPIID_AC( string strfilePPIMapp )
{

	ofstream outf13("parse_ProteinIPIID_AC2008.txt");

	//	CMassPeptideIons* pMassPeptideIons=NULL;
	ifstream is (strfilePPIMapp.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strfilePPIMapp).c_str());
	}
	int irow=0;
	int iPeptideRow=0;


	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{


		istringstream strStrem(buffer);

		string strProtienIPIID;
		strStrem>>strProtienIPIID;

		string strProteinACID;
		strStrem>>strProteinACID;


		string::size_type pos2 = strProteinACID.find(';');
		string strSwissProtIDRest=strProteinACID;
		while (pos2!=-1)
		{
			string strSingleProtein=strSwissProtIDRest.substr(0, pos2);
			outf13<<strProtienIPIID<<'\t'<<strSingleProtein<<endl;


			strSwissProtIDRest=strSwissProtIDRest.substr(pos2+1, strSwissProtIDRest.length());
			pos2 = strSwissProtIDRest.find(';');

		}


		outf13<<strProtienIPIID<<'\t'<<strSwissProtIDRest<<endl;


	}
	cout<<irow;
	cin>>irow;
}
