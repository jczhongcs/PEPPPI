#include "stdafx.h"
#include "Protein.h"

#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <set>
#include <xfunctional>
#include <algorithm>

string&   replace_all(string&   str,const   string&   old_value,const   string&   new_value)   
{   
	while(true)   {   
		string::size_type   pos(0);   
		if(   (pos=str.find(old_value))!=string::npos   )   
			str.replace(pos,old_value.length(),new_value);   
		else   break;   
	}   
	return   str;   
}   


CProtein::CProtein(string& strHead,string& strSequence)
{
	//parse the head and set sequence
	//
	m_ivalidTheoryPeptideCount=0;
	m_strHead=strHead;
	string::size_type pos2 = strHead.find(':');
	string strProteinNameRest=replace_all(strHead," ","|");
	while (pos2!=-1)
	{
		string strProteinTitle=strProteinNameRest.substr(0, pos2);
//		outf13<<strSingleProtein<<'\t'<<strSpectrumName<<'\t'<<dXcorScore<<'\t'<<strPeptideSequence<<'\t'<<dPeptideMass<<'\t'<<iCharge<<endl;
		string::size_type pos3 = strProteinNameRest.find('|');
		string strNames=strProteinNameRest.substr(pos2+1, pos3-pos2-1);
		if (strProteinTitle==">IPI")
		{
			string::size_type pos3 = strProteinNameRest.find('.');
			
			m_strIPIName=strProteinNameRest.substr(pos2+1, pos3-pos2-1);

			m_strProteinName=m_strIPIName;

		}else if (strProteinTitle=="SWISS-PROT")
		{
			ExtractNames(strNames, m_strSwissName);

		}else if (strProteinTitle=="TREMBL")
		{
			ExtractNames(strNames, m_strTREMBL);

		}else if (strProteinTitle=="ENSEMBL")
		{
			ExtractNames(strNames, m_strENSEMBL);

		}else if (strProteinTitle=="REFSEQ")
		{
			ExtractNames(strNames, m_strREFSEQ);

		}

		pos2 = strProteinNameRest.find('|');
		//////////////////////////////////////////////////////////////////////////
		//no cut
		if (pos2==-1)
		{

			break;
			
		}
		strProteinNameRest=strProteinNameRest.substr(pos2+1, strProteinNameRest.length());
		pos2 = strProteinNameRest.find(':');

	}
	//////////////////////////////////////////////////////////////////////////
	for (int i=0;i<m_strSwissName.size();i++)
	{
		string::size_type pos = m_strSwissName[i].find('-');
		
		if (pos!=-1)
		{
			m_strSwissName[i]=m_strSwissName[i].substr(0, pos);
		}
	}


	m_strSequence=strSequence;
	string strEn("Trypsin - N-term K or R, not C-term P");
	objEnzym.SplitProteinWithEnzyme(m_strSequence,strEn,m_vTheoryPeptide,m_ivalidTheoryPeptideCount);
	m_dProbability=0.0;
	m_dSupportValueFromOtherProtein=0.0;

}
//////////////////////////////////////////////////////////////////////////
//处理yeast序列数据
CProtein::CProtein( string& strHead,string& strSequence,bool bYeast )
{
	m_ivalidTheoryPeptideCount=0;
	istringstream strStrem(strHead);
	string strYeastName,strGeneName;
	strStrem>>strYeastName>>strGeneName;
	m_strProteinName=m_strIPIName=strYeastName.substr(1,strYeastName.length());//此处为与人类兼容
	m_strSequence=strSequence;
	m_strHead=strHead;
	m_strGeneName=strGeneName;

	string strEn("Trypsin - N-term K or R, not C-term P");
	objEnzym.SplitProteinWithEnzyme(m_strSequence,strEn,m_vTheoryPeptide,m_ivalidTheoryPeptideCount);
	m_dProbability=0.0;
	m_dSupportValueFromOtherProtein=0.0;
}
//////////////////////////////////////////////////////////////////////////
//处理human sequence data
CProtein::CProtein(string& strHead,string& strSequence,bool bHuman,bool bWithDecoy)
{
	//parse the head and set sequence
	//
	m_ivalidTheoryPeptideCount=0;
	m_strHead=strHead;

	string::size_type pos2 = strHead.find(">IPI");
	if (pos2!=-1)//target
	{
	
		//replace all space to |
		string strProteinNameRest=replace_all(strHead," ","|");
		//>IPI00000001 IPI:IPI00000001.2|SWISS-PROT:O95793-1|TREMBL:A8K622;Q59F99|ENSEMBL:ENSP00000360922;ENSP00000379466|REFSEQ:NP_059347|H-INV:HIT000329496|VEGA:OTTHUMP00000031233 Tax_Id=9606 Gene_Symbol=STAU1 Isoform Long of Double-stranded RNA-binding protein Staufen homolog 1
		//remove the >IPI00000001 
		pos2 = strProteinNameRest.find('|');
		strProteinNameRest=strProteinNameRest.substr(pos2+1, strProteinNameRest.length());
		pos2 = strProteinNameRest.find(':');

		while (pos2!=-1)
		{

			string strProteinTitle=strProteinNameRest.substr(0, pos2);
			//		outf13<<strSingleProtein<<'\t'<<strSpectrumName<<'\t'<<dXcorScore<<'\t'<<strPeptideSequence<<'\t'<<dPeptideMass<<'\t'<<iCharge<<endl;
			string::size_type pos3 = strProteinNameRest.find('|');
			string strNames=strProteinNameRest.substr(pos2+1, pos3-pos2-1);
			if (strProteinTitle=="IPI")
			{
				string::size_type pos3 = strProteinNameRest.find('.');

				m_strIPIName=strProteinNameRest.substr(pos2+1, pos3-pos2-1);

				m_strProteinName=m_strIPIName;

			}else if (strProteinTitle=="SWISS-PROT")
			{
				ExtractNames(strNames, m_strSwissName);

			}else if (strProteinTitle=="TREMBL")
			{
				ExtractNames(strNames, m_strTREMBL);

			}else if (strProteinTitle=="ENSEMBL")
			{
				ExtractNames(strNames, m_strENSEMBL);

			}else if (strProteinTitle=="REFSEQ")
			{
				ExtractNames(strNames, m_strREFSEQ);

			}

			pos2 = strProteinNameRest.find('|');
			//////////////////////////////////////////////////////////////////////////
			//no cut
			if (pos2==-1)
			{

				break;

			}
			strProteinNameRest=strProteinNameRest.substr(pos2+1, strProteinNameRest.length());
			pos2 = strProteinNameRest.find(':');

		}
/*
		//////////////////////////////////////////////////////////////////////////
		//remove swissPort name - number
		for (int i=0;i<m_strSwissName.size();i++)
		{
			string::size_type pos = m_strSwissName[i].find('-');

			if (pos!=-1)
			{
				m_strSwissName[i]=m_strSwissName[i].substr(0, pos);
			}
		}
*/
		m_bDecoy=false;
	}
	pos2 = strHead.find(">DECOY_");
	if (pos2!=-1)//Decoy
	{
		string strProteinNameRest=replace_all(strHead," ","|");
		string::size_type pos1 = strProteinNameRest.find('>');
		pos2 = strProteinNameRest.find('|');
		m_strIPIName=strProteinNameRest.substr(pos1+1, pos2-pos1-1);
		m_strProteinName=m_strIPIName;
		m_bDecoy=true;

	}

	m_strSequence=strSequence;
	string strEn("Trypsin - N-term K or R, not C-term P");
	objEnzym.SplitProteinWithEnzyme(m_strSequence,strEn,m_vTheoryPeptide,m_ivalidTheoryPeptideCount);
	m_dProbability=0.0;
	m_dSupportValueFromOtherProtein=0.0;

}

CProtein::CProtein( string& strHead,string& strSequence,bool bHuman,bool bWithDecoy,bool ensembl )
{
	//parse the head and set sequence
	//
	m_ivalidTheoryPeptideCount=0;
	m_strHead=strHead;
	string::size_type pos2 = strHead.find('>');
	string strProteinNameRest=replace_all(strHead," ","|");
	string::size_type pos3 = strProteinNameRest.find('|');
	m_strIPIName=strProteinNameRest.substr(pos2+1, pos3-pos2-1);
	m_strProteinName=m_strIPIName;
	m_strENSEMBL.push_back(m_strProteinName);
	m_strSequence=strSequence;
	string strEn("Trypsin - N-term K or R, not C-term P");
	objEnzym.SplitProteinWithEnzyme(m_strSequence,strEn,m_vTheoryPeptide,m_ivalidTheoryPeptideCount);
	m_dProbability=0.0;
	m_dSupportValueFromOtherProtein=0.0;
}

CProtein::~CProtein(void)
{
}

void CProtein::parse_Sequence( CPeptideProtein& peptideProtein )
{
//////////////////////////////////////////////////////////////////////////
//int type, mapPeptideProtein

}

void CProtein::ExtractNames( string &strName, vector<string>& vstr )
{
	string::size_type pos = strName.find(';');
	string strNameRest=strName;
	while (pos!=-1)
	{
		string strSingleName=strNameRest.substr(0, pos);
		vstr.push_back(strSingleName);


		strNameRest=strNameRest.substr(pos+1, strNameRest.length());
		pos = strNameRest.find(';');

	}	
	vstr.push_back(strNameRest);
}

bool CProtein::findProtein( string str )
{
/*
	if (m_strHead.find(str)!=-1)
	{
		return true;
	}
*/
	if (m_strIPIName==str) return true;
	if(find(m_strSwissName.begin(),m_strSwissName.end(),str)!=m_strSwissName.end()) return true;
	if(find(m_strTREMBL.begin(),m_strTREMBL.end(),str)!=m_strTREMBL.end()) return true;
	if(find(m_strENSEMBL.begin(),m_strENSEMBL.end(),str)!=m_strENSEMBL.end()) return true;
	if(find(m_strREFSEQ.begin(),m_strREFSEQ.end(),str)!=m_strREFSEQ.end()) return true;



	return false;
}

void CProtein::outFileProteinName( ofstream& out )
{
	int isize=m_strSwissName.size()+m_strTREMBL.size()+m_strENSEMBL.size()+m_strREFSEQ.size();
	if (isize==0)
	{
		out<<m_strIPIName<<"\tIPINAME\t"<<m_strIPIName<<endl;
		return;

	}
	for(int i=0;i<m_strSwissName.size();i++)
		out<<m_strIPIName<<"\tSwissID\t"<<m_strSwissName[i]<<endl;
	for(int i=0;i<m_strTREMBL.size();i++)
		out<<m_strIPIName<<"\tTREMBLID\t"<<m_strTREMBL[i]<<endl;
	for(int i=0;i<m_strENSEMBL.size();i++)
		out<<m_strIPIName<<"\tENSEMBLID\t"<<m_strENSEMBL[i]<<endl;
	for(int i=0;i<m_strREFSEQ.size();i++)
		out<<m_strIPIName<<"\tREFSEQID\t"<<m_strREFSEQ[i]<<endl;
	

}

void CProtein::outFilePeptideProbability( ofstream& out,double myValue  ,bool bShowDetail)
{
//	find_if(m_vSupportPeptide.begin(),m_vSupportPeptide.end(),
	
//	double myValue=0.99;
	int iConfidence = count_if(m_vSupportPeptide.begin(), m_vSupportPeptide.end(),[myValue](CMassPeptideXcorr*  x) { return x->m_dProbability > myValue; }); 
	out<< iConfidence<<'\t';
//	int dayu5=count_if(m_vSupportPeptide.begin(),m_vSupportPeptide.end(),bind2nd(greater<double>(),5.09));
	if (bShowDetail)
		for (CMassPeptideXcorr* pMP:m_vSupportPeptide)
		{
			out<<pMP->m_dProbability/*<<'\t'*/<<(pMP->m_bUnique?'U':'S')<<'\t';
		}
	
}

void CProtein::calculateProbabilitySharedOne()
{
	double dProbability=1.0;
	for (CMassPeptideXcorr* pPeptideXcorr:m_vSupportPeptide)
	{
		dProbability=dProbability*(1-pPeptideXcorr->m_dProbability);
	}
	m_dProbability=pow((1-dProbability),-log(m_vSupportPeptide.size()/m_ivalidTheoryPeptideCount));
}

void CProtein::calculateProbabilitySharedDivide()
{
	double dProbability=1.0;
	for (CMassPeptideXcorr* pPeptideXcorr:m_vSupportPeptide)
	{
		int iPeptideSharedProtein=1;
		if (pPeptideXcorr->m_bUnique==false)//shared peptide 均分
		{
			iPeptideSharedProtein=pPeptideXcorr->m_vProtein.size();
		}

		dProbability=dProbability*(1-pPeptideXcorr->m_dProbability/iPeptideSharedProtein);
		
	}
	m_dProbability=pow((1-dProbability),-log(m_vSupportPeptide.size()/m_ivalidTheoryPeptideCount));
}

/*

void CProtein::calculateProbabilitySharedWeighted(CListMassPeptideXcorr* pListMassPeptideXcorr,CListProtein* pListProtein)
{
	double dProbability=1.0;
	for (CMassPeptideXcorr* pPeptideXcorr:m_vSupportPeptide)//Loop SupportPeptide
	{
		int iPeptideSharedProtein=1;
		double dOriginPeptideSP=1.0;
		if (pPeptideXcorr->m_bUnique==false)//shared peptide 
		{
			double dOtherProteinPeptideValue=0.0;
			dOriginPeptideSP=getSameMassSPValue(pListMassPeptideXcorr, pPeptideXcorr);//get the peptide sigma sp value, mass match to peptide 

			vector<double> vdPeptideSP(pPeptideXcorr->m_vProtein.size(),0.0);//save the sibling peptide value
			for (int j=0;j<pPeptideXcorr->m_vProtein.size();j++)//Loop share peptide to Protein
			{
				vector<CMassPeptideXcorr*> vMassPeptide = pListProtein->getSupportPeptideFromProteinStr(pPeptideXcorr->m_vProtein[j]);
				double dPeptideSP=0.0;
				for (int j1=0;j1<vMassPeptide.size();j1++)//every protein has some peptides 
				{
					if (vMassPeptide[j1]->m_strPeptideSequence!=pPeptideXcorr->m_strPeptideSequence)//except for the original peptideSequence
					{
						dPeptideSP += getSameMassSPValue(pListMassPeptideXcorr, vMassPeptide[j1]);
					}

				}
				dPeptideSP/=(pPeptideXcorr->m_vProtein.size()-1);

				vdPeptideSP.push_back(dPeptideSP);
				if (pPeptideXcorr->m_vProtein[j]!=m_strIPIName)
				{
					dOtherProteinPeptideValue+=dPeptideSP;
				}

			}
			dOriginPeptideSP=abs(dOriginPeptideSP-dOtherProteinPeptideValue)/dOriginPeptideSP;


		}

		dProbability=dProbability*(1-dOriginPeptideSP);//no exist probability

	}
	dProbability=pow((1-dProbability),-log(m_vSupportPeptide.size()/m_vTheoryPeptide.size()));
}
*/

void CProtein::calculateProbabilitySharedLogisticRegression()
{
/*
	double dProbability=1.0;
	for (CMassPeptideXcorr* pPeptideXcorr:m_vSupportPeptide)//Loop SupportPeptide
	{
		int iPeptideSharedProtein=1;
		double dOriginPeptideSP=1.0;
		if (pPeptideXcorr->m_bUnique==false)//shared peptide 
		{
			double dOtherProteinPeptideValue=0.0;
			dOriginPeptideSP=getSameMassSPValue(pListMassPeptideXcorr, pPeptideXcorr);//get the peptide sigma sp value, mass match to peptide 

			vector<double> vdPeptideSP(pPeptideXcorr->m_vProtein.size(),0.0);//save the sibling peptide value
			for (int j=0;j<pPeptideXcorr->m_vProtein.size();j++)//Loop share peptide to Protein
			{
				vector<CMassPeptideXcorr*> vMassPeptide = pListProtein->getSupportPeptideFromProteinStr(pPeptideXcorr->m_vProtein[j]);
				double dPeptideSP=0.0;
				for (int j1=0;j1<vMassPeptide.size();j1++)//every protein has some peptides 
				{
					if (vMassPeptide[j1]->m_strPeptideSequence!=pPeptideXcorr->m_strPeptideSequence)//except for the original peptideSequence
					{
						dPeptideSP += getSameMassSPValue(pListMassPeptideXcorr, vMassPeptide[j1]);
					}

				}
				dPeptideSP/=(pPeptideXcorr->m_vProtein.size()-1);

				vdPeptideSP.push_back(dPeptideSP);
				if (pPeptideXcorr->m_vProtein[j]!=m_strIPIName)
				{
					dOtherProteinPeptideValue+=dPeptideSP;
				}

			}
			dOriginPeptideSP=abs(dOriginPeptideSP-dOtherProteinPeptideValue)/dOriginPeptideSP;


		}

		dProbability=dProbability*(1-dOriginPeptideSP);//no exist probability

	}
	dProbability=pow((1-dProbability),-log(m_vSupportPeptide.size()/m_vTheoryPeptide.size()));*/
}

double CProtein::getSameMassSPValue( CListMassPeptideXcorr* pListMassPeptideXcorr, CMassPeptideXcorr* pPeptideXcorr )
{
	double dPeptideSP=0.0;
	for (int i=0;i<pListMassPeptideXcorr->m_mapvSameMassPeptides[pPeptideXcorr->m_strPeptideSequence].size();i++)
	{
		dPeptideSP+=pListMassPeptideXcorr->m_vpMassPeptideXcorr[pListMassPeptideXcorr->m_mapvSameMassPeptides[pPeptideXcorr->m_strPeptideSequence][i]]->m_dSp;
	}	return dPeptideSP;
}

void CListProtein::readparse_Yeastfile( string strFileName )
{
	//	CMassPeptideIons* pMassPeptideIons=NULL;


	//outf13 can be out to file


	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}


	int iProteinCount=0;
	string strProteinHead="";
	string strProteinSequence="";
	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{

		if (buffer[0]=='>')
		{
			if (strProteinHead.length()>1)
			{
				CProtein* objProtein=new CProtein(strProteinHead,strProteinSequence,true);
				vProteinList.push_back(objProtein);
				mapListPos.insert(make_pair(objProtein->m_strIPIName,iProteinCount));
				//				outf13<<"\t"<<objProtein.m_strIPIName<<"\t"<<objProtein.m_strREFSEQ
				iProteinCount++;
			}

			//////////////////////////////////////////////////////////////////////////
			//add protein



			//////////////////////////////////////////////////////////////////////////
			//new proteinHead
			strProteinHead=buffer;
			strProteinSequence="";

			continue;
		}
		strProteinSequence+=buffer;


	}
	CProtein* objProtein=new CProtein(strProteinHead,strProteinSequence,true);
	vProteinList.push_back(objProtein);
	mapListPos.insert(make_pair(objProtein->m_strIPIName,iProteinCount));
	iProteinCount++;


	outProteinAndSequence(objProtein);


}
void CListProtein::readparse_Esembelfile( string strFileName )
{

	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}


	int iProteinCount=0;
	string strProteinHead="";
	string strProteinSequence="";
	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{

		if (buffer[0]=='>')
		{
			if (strProteinHead.length()>1)
			{
				CProtein* objProtein=new CProtein(strProteinHead,strProteinSequence,true,true,true);
				vProteinList.push_back(objProtein);
				mapListPos.insert(make_pair(objProtein->m_strIPIName,iProteinCount));//这里的IPI实质是Esembl
				//				outf13<<"\t"<<objProtein.m_strIPIName<<"\t"<<objProtein.m_strREFSEQ
				iProteinCount++;
			}

			//////////////////////////////////////////////////////////////////////////
			//add protein



			//////////////////////////////////////////////////////////////////////////
			//new proteinHead
			strProteinHead=buffer;
			strProteinSequence="";

			continue;
		}
		strProteinSequence+=buffer;


	}
	CProtein* objProtein=new CProtein(strProteinHead,strProteinSequence,true,true);
	vProteinList.push_back(objProtein);
	mapListPos.insert(make_pair(objProtein->m_strIPIName,iProteinCount));
	iProteinCount++;
}

void CListProtein::readparse_file( string strFileName )
{
	//	CMassPeptideIons* pMassPeptideIons=NULL;


	//outf13 can be out to file


	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}

	
	int iProteinCount=0;
	string strProteinHead="";
	string strProteinSequence="";
	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{

		if (buffer[0]=='>')
		{
			if (strProteinHead.length()>1)
			{
				CProtein* objProtein=new CProtein(strProteinHead,strProteinSequence,true,true);
				vProteinList.push_back(objProtein);
				mapListPos.insert(make_pair(objProtein->m_strIPIName,iProteinCount));
//				outf13<<"\t"<<objProtein.m_strIPIName<<"\t"<<objProtein.m_strREFSEQ
				iProteinCount++;
			}

			//////////////////////////////////////////////////////////////////////////
			//add protein



			//////////////////////////////////////////////////////////////////////////
			//new proteinHead
			strProteinHead=buffer;
			strProteinSequence="";

			continue;
		}
		strProteinSequence+=buffer;


	}
	CProtein* objProtein=new CProtein(strProteinHead,strProteinSequence,true,true);
	vProteinList.push_back(objProtein);
	mapListPos.insert(make_pair(objProtein->m_strIPIName,iProteinCount));
	iProteinCount++;

}

ifstream& CListProtein::getcleanline( ifstream& is, string& buffer )
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

CListProtein::~CListProtein( void )
{
	for (int i =0;i<vProteinList.size();i++)
	{
		delete vProteinList[i];
	}
	vProteinList.clear();
}

void CListProtein::OutFileProteinListHeadInfo( string strFileName )
{

	ofstream out13("ProteinListInfo");
	for (int i =0;i<vProteinList.size();i++)
	{
		vProteinList[i]->outFileProteinName(out13);
	}

}

void CListProtein::MatchProteinHeadList( string str,vector<string>& vIPIName )
{
	for (int i =0;i<vProteinList.size();i++)
	{

		if (vProteinList[i]->findProtein(str))
		{
			vIPIName.push_back(vProteinList[i]->m_strIPIName);
		}
	}

}


void CListProtein::mergePPI( string str )
{
	ofstream out("mergePPI");
	ifstream is (str.c_str());
	if(!is)
	{
		throw string(("Can't open file"+str).c_str());
	}
	string s1;
	string s2;
	string s3;
	string s4;
	string s5;
	vector<string> vString;
	while(is>>s1)
	{
		is>>s2;
		is>>s3;
		s4=s1+s2;
		s5=s2+s1;
		if ((find(vString.begin(),vString.end(),s4)==vString.end())&&(find(vString.begin(),vString.end(),s5)==vString.end()))
		{
			vString.push_back(s4);
			out<<s1<<'\t'<<s2<<'\t'<<s3<<endl;
		}
		
	}

}

void CListProtein::outProteinAndSequence( CProtein* objProtein ,int iCount/*=100*/ )
{
	ofstream* out13=NULL;
	for (int i =0;i<vProteinList.size();i++)
	{
		if(i%iCount==0)
		{
			if (out13)
			{
				out13->flush();
				out13->close();
				delete out13;

			}
			stringstream ss;
			ss<<i/iCount;
			out13=new ofstream("protein"+ss.str());
		}
		(*out13) << '>'<< vProteinList[i]->m_strIPIName<<'\n'<<vProteinList[i]->m_strSequence<<endl;

	}
}

//////////////////////////////////////////////////////////////////////////
//读入参考集
void CListProtein::readGoldenSet( string strFileName )
{

	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}
	string s1;

	while(is>>s1)
	{
	
		vYeastGoldenSet.push_back(s1);
	}

}
//////////////////////////////////////////////////////////////////////////
//查找是不是在参考集中
int CListProtein::isInGoldenset( string proteinName )
{
	if ( find(vYeastGoldenSet.begin(),vYeastGoldenSet.end(),proteinName)!=vYeastGoldenSet.end())
	{
		return 1;
	}
	return 0;
}

bool CListProtein::isUniquePeptide( CMassPeptideXcorr* pPeptide )
{
	return false;
}


bool CListProtein::OutNoRepeatPeptideInProtein( CMassPeptideXcorr* pPeptide,ofstream& out )
{
	set<string> setProtein;
	for (int j = 0;j< pPeptide->m_vProtein.size();j++)
	{
		int i=iGetProteinPos( pPeptide->m_vProtein[j]);
		if (i>-1)//finded
			setProtein.insert(vProteinList[i]->m_strIPIName);
	
		
	}


	for (string strProName:setProtein)
	{
		out << pPeptide->m_strPeptideSequence<<'\t'<<strProName<<'\t'<<pPeptide->m_dProbability<<"\tFindS"<<endl;
	
	}
	return true;
}

bool CListProtein::isPeptideInProtein( CMassPeptideXcorr* pPeptide,ofstream& out )
{
	set<int> setProtein;
	
	for (int j = 0;j< pPeptide->m_vProtein.size();j++)
	{
		int i=iGetProteinPos( pPeptide->m_vProtein[j]);
		if (i>-1)//finded
		{

			if (vProteinList[i]->findProteinSequence(pPeptide->m_strPeptideSequence))
			{ 
				out << pPeptide->m_strPeptideSequence<<'\t'<<vProteinList[i]->m_strIPIName<<'\t'<<pPeptide->m_dProbability<<"\tFindS"<<endl;
				setProtein.insert(i);
			}
			else
			{
				out << pPeptide->m_strPeptideSequence<<'\t'<<vProteinList[i]->m_strIPIName<<'\t'<<pPeptide->m_dProbability<<"\tNoFindS"<<endl;
//				return false;
			}
			
		}
		else
		{
			out << pPeptide->m_strPeptideSequence<<'\t'<<pPeptide->m_vProtein[j]<<'\t'<<pPeptide->m_dProbability<<"\tNoMacthProtein"<<endl;
//			return false;
		}
	}
	return true;
}

int CListProtein::iGetProteinPos( string proteinName )
{
	bool bDecoy=false;
	if (proteinName.find("DECOY_")!=-1)
	{
		bDecoy=true;
		proteinName=replace_all(proteinName,"DECOY_","");
	}

	for (int i =0;i<vProteinList.size();i++)
	{
		if (vProteinList[i]->m_strHead.find(proteinName)!=-1)
		{
			if (bDecoy)
			{
				return i+1;//return decoy
			}
			return i;
		}
// 
// 		if (vProteinList[i]->findProtein(proteinName))
// 		{
// 			return i;
// 		}
	}
	return -1;
}

bool CProtein::findProteinSequence( string str )
{

	if (m_strSequence.find(str)!=-1)
	{
		return true;
	}
}