#include "stdafx.h"
#include "ListMassPeptideXcorr.h"


#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
typedef pair<string, vector<int> > PROTEINPAIR;


#define MassProbability 0.99;


int cmp(const PROTEINPAIR& x,const PROTEINPAIR& y)
{
	return x.second[0]>y.second[0];
}


bool getHighConfidencePeptides(CMassPeptideXcorr* pMassPeptideXcorr)
{
	return pMassPeptideXcorr->m_dProbability > MassProbability;
}


CListMassPeptideXcorr::CListMassPeptideXcorr(void)
{
}


CListMassPeptideXcorr::~CListMassPeptideXcorr(void)
{
	for (vector< CMassPeptideXcorr* >::iterator it=m_vpMassPeptideXcorr.begin();
		it!=m_vpMassPeptideXcorr.end();
		it++)
	{
		delete *it;
	}
	m_vpMassPeptideXcorr.clear();
}
//////////////////////////////////////////////////////////////////////////
//read from the file
void CListMassPeptideXcorr::readparse_file( string strFileName ,int iType)
{
	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}
	int irow=0;

	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
		string strProteinName;
		string strSpectrum;
		double dDeltaCN;
		string strPeptideSequence;
		double dMZ;
		double dCalculateH;
		double dActuralPeptide;
		int iSpectrumCharge;
		bool bUnique=true;
		string strFileName;
		int iScanNO;
		double dRetenTime;
		int iMissCleave;
		double dMassDiff;
		double dProbability;
		double dSp;

		bool bIPIDecoy;
		char cPreResidue;
		char cPostResidue;
		string strFullMod;
		int iMatch;
		int iPredict;
		double dXCorr=0.0;

		istringstream strStrem(buffer);//×Ö·û´®±ä³ÉÁ÷

				strStrem>>strProteinName>>strSpectrum>>dXCorr>>dDeltaCN>>strPeptideSequence>>dMZ>>dCalculateH>>dActuralPeptide>>iSpectrumCharge;
//		strStrem>>strProteinName>>strSpectrum>>dXCorr>>dDeltaCN>>strPeptideSequence>>dMZ>>dCalculateH>>dActuralPeptide>>iSpectrumCharge;
		//Peptide corresponding proteins
		string::size_type pos2 = strProteinName.find(',');
		if (pos2!=-1)	bUnique=false;

		string strProteinNameRest=strProteinName;
		while (pos2!=-1)
		{
			string strSingleProtein=strProteinNameRest.substr(0, pos2);

			CMassPeptideXcorr* pMassPeptideXcorr=new CMassPeptideXcorr(strSingleProtein,strSpectrum,dXCorr,dDeltaCN,strPeptideSequence,dMZ,dCalculateH,dActuralPeptide,iSpectrumCharge,bUnique);
			//Protein accession numbers	Spectrum name	SEQUEST XCorr score	Peptide sequence	Observed m/z	Calculated +1H Peptide Mass (AMU)	Actual peptide mass (AMU)	Spectrum charge

			m_vpMassPeptideXcorr.push_back(pMassPeptideXcorr);


			strProteinNameRest=strProteinNameRest.substr(pos2+1, strProteinNameRest.length());
			pos2 = strProteinNameRest.find(',');

		}
		CMassPeptideXcorr* pMassPeptideXcorr=new CMassPeptideXcorr(strProteinNameRest,strSpectrum,dXCorr,dDeltaCN,strPeptideSequence,dMZ,dCalculateH,dActuralPeptide,iSpectrumCharge,bUnique);
		//Protein accession numbers	Spectrum name	SEQUEST XCorr score	Peptide sequence	Observed m/z	Calculated +1H Peptide Mass (AMU)	Actual peptide mass (AMU)	Spectrum charge

		m_vpMassPeptideXcorr.push_back(pMassPeptideXcorr);


	}
}

void CListMassPeptideXcorr::readparseFidoe_file( string strFileName )
{
	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}
	int irow=0;
	ofstream out("fidoResultMaxdata02gt05a1.txt");
	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
		string strProteinName;
		string strProteinName1;
		string strProteinName2;
// 		string strProteinName3;
// 		string strProteinName4;
// 		string strProteinName5;
		
		double dProbability;
		istringstream strStrem(buffer);
		strStrem>>dProbability>>strProteinName>>strProteinName1>>strProteinName2/*>>strProteinName3>>strProteinName4>>strProteinName5*/;
		out<<strProteinName<<'\t'<<dProbability<<endl;
		out<<strProteinName1<<'\t'<<dProbability<<endl;
		out<<strProteinName2<<'\t'<<dProbability<<endl;
// 		out<<strProteinName3<<'\t'<<dProbability<<endl;
// 		out<<strProteinName4<<'\t'<<dProbability<<endl;
// 		out<<strProteinName5<<'\t'<<dProbability<<endl;
// 		

	}
}

void CListMassPeptideXcorr::readparseFidoe_file( string strFileName,char charSplitChar )
{
	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}
	int irow=0;
	ofstream out("fidoResultMaxdata02gt05a1.txt");
	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
		string strProteinName;
		string strProteinName1;
		string strProteinName2;
// 		string strProteinName3;
// 		string strProteinName4;
// 		string strProteinName5;
		
		double dProbability;
		istringstream strStrem(buffer);
		strStrem>>dProbability>>strProteinName;
		string::size_type pos2 = strProteinName.find(charSplitChar);
		//		if (pos2!=-1)	bUnique=false;

		string strProteinNameRest=strProteinName;
		while (pos2!=-1)
		{
			string strSingleProtein=strProteinNameRest.substr(0, pos2);


			out<<strSingleProtein<<'\t'<<dProbability<<endl;


			strProteinNameRest=strProteinNameRest.substr(pos2+1, strProteinNameRest.length());
			pos2 = strProteinNameRest.find(charSplitChar);

		}

		out<<strProteinNameRest<<'\t'<<dProbability<<endl;

	}
}

//////////////////////////////////////////////////////////////////////////
//read from the file
void CListMassPeptideXcorr::readparsePrePareSimple_file( string strFileName,char charSplitChar )
{
	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}
	int irow=0;

	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
		string strProteinName;
		string strSpectrum;
		double dProbability;
		istringstream strStrem(buffer);
		strStrem>>strSpectrum>>strProteinName>>dProbability;
		string strPepIdentify;
		stringstream ss;
		ss<<dProbability;

		strPepIdentify=strSpectrum+ss.str();

		string::size_type pos2 = strProteinName.find(charSplitChar);
//		if (pos2!=-1)	bUnique=false;

		string strProteinNameRest=strProteinName;
		while (pos2!=-1)
		{
			string strSingleProtein=strProteinNameRest.substr(0, pos2);

			

			if(m_mapMassPeptidePos.find(strPepIdentify)==m_mapMassPeptidePos.end())
			{
				CMassPeptideXcorr* pMassPeptideXcorr=new CMassPeptideXcorr( strSingleProtein, strSpectrum, dProbability);

				pMassPeptideXcorr->m_vProtein.push_back(strSingleProtein);
				pMassPeptideXcorr->m_mapProteinPos.insert(make_pair(strSingleProtein,pMassPeptideXcorr->m_vProtein.size()-1));

				m_vpMassPeptideXcorr.push_back(pMassPeptideXcorr);

				m_mapMassPeptidePos.insert(make_pair(strPepIdentify,irow));

				//////////////////////////////////////////////////////////////////////////
				//insert the same peptidesequence from difference mass data
				if( m_mapvSameMassPeptides.find(strSpectrum)==m_mapvSameMassPeptides.end())
				{
					vector<int> v;
					v.push_back(irow);
					m_mapvSameMassPeptides.insert(make_pair(strSpectrum,v));
				}
				else
				{
					m_mapvSameMassPeptides[strSpectrum].push_back(irow);
				}
				irow++;



			}
			else
			{
				if(find(m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.begin(),m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.end(),strSingleProtein)==m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.end())
				{
					m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.push_back(strSingleProtein);
					m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_mapProteinPos.insert(make_pair(strSingleProtein,m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.size()-1));
				}
			}

//			m_vpMassPeptideXcorr.push_back(pMassPeptideXcorr);


			strProteinNameRest=strProteinNameRest.substr(pos2+1, strProteinNameRest.length());
			pos2 = strProteinNameRest.find(charSplitChar);

		}
		//add the rest protein
		if(m_mapMassPeptidePos.find(strPepIdentify)==m_mapMassPeptidePos.end())
		{
			CMassPeptideXcorr* pMassPeptideXcorr=new CMassPeptideXcorr( strProteinNameRest, strSpectrum, dProbability);

			pMassPeptideXcorr->m_vProtein.push_back(strProteinNameRest);
			pMassPeptideXcorr->m_mapProteinPos.insert(make_pair(strProteinNameRest,pMassPeptideXcorr->m_vProtein.size()-1));

			m_vpMassPeptideXcorr.push_back(pMassPeptideXcorr);

			m_mapMassPeptidePos.insert(make_pair(strPepIdentify,irow));

			//////////////////////////////////////////////////////////////////////////
			//insert the same peptidesequence from difference mass data
			if( m_mapvSameMassPeptides.find(strSpectrum)==m_mapvSameMassPeptides.end())
			{
				vector<int> v;
				v.push_back(irow);
				m_mapvSameMassPeptides.insert(make_pair(strSpectrum,v));
			}
			else
			{
				m_mapvSameMassPeptides[strSpectrum].push_back(irow);
			}
			irow++;



		}
		else
		{
			if(find(m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.begin(),m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.end(),strProteinNameRest)==m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.end())
			{
				m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.push_back(strProteinNameRest);
				m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_mapProteinPos.insert(make_pair(strProteinNameRest,m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.size()-1));
			}
		}
	}
	//set average value to protein which have some shared peptide

	//////////////////////////////////////////////////////////////////////////
	//fidoFormatOutput
	ofstream outfido("fidoOutSequestData2MaxGt05.txt");
	for (int i=0;i<m_vpMassPeptideXcorr.size();i++)
	{
		outfido<<'e'<<' '<<m_vpMassPeptideXcorr[i]->m_strSpectrum<<endl;
		int iPeptideProteinSize=m_vpMassPeptideXcorr[i]->m_vProtein.size();
		for (int j = 0;j< iPeptideProteinSize;j++)
		{
			outfido<<'r'<<' '<<m_vpMassPeptideXcorr[i]->m_vProtein[j]<<endl;

		}
		outfido<<'p'<<' '<< m_vpMassPeptideXcorr[i]->m_dProbability<<endl;
	}
	//////////////////////////////////////////////////////////////////////////
	//peptideProteinOutput
	ofstream outpeptideProteinOutput("peptideProteinOutputData2.txt");
	for (int i=0;i<m_vpMassPeptideXcorr.size();i++)
	{
		int iPeptideProteinSize=m_vpMassPeptideXcorr[i]->m_vProtein.size();
		for (int j = 0;j< iPeptideProteinSize;j++)
		{
			outpeptideProteinOutput<<m_vpMassPeptideXcorr[i]->m_strSpectrum<<'\t'<<m_vpMassPeptideXcorr[i]->m_vProtein[j]<<'\t'<< m_vpMassPeptideXcorr[i]->m_dProbability<<endl;

		}
	}
	SetProbPeptideToProteinOnAverage();
	//SetProbPeptideToProteinOnUniqueNoConsiderShared();
}

//////////////////////////////////////////////////////////////////////////
//read from the file
void CListMassPeptideXcorr::convertToFido( string strFileName )
{
	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}
	int irow=0;

	ofstream outfido("fidoOutSequestdata02gt05.txt");

	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
		string strProteinName;
		string strSpectrum;
		double dProbability;
		istringstream strStrem(buffer);
		strStrem>>strSpectrum>>strProteinName>>dProbability;
		string strPepIdentify;
		stringstream ss;
		ss<<dProbability;

		string::size_type pos2 = strProteinName.find(',');



		outfido<<'e'<<' '<<strSpectrum<<endl;
		
		

		string strProteinNameRest=strProteinName;
		while (pos2!=-1)
		{
			string strSingleProtein=strProteinNameRest.substr(0, pos2);


			outfido<<'r'<<' '<<strSingleProtein<<endl;


			

			strProteinNameRest=strProteinNameRest.substr(pos2+1, strProteinNameRest.length());
			pos2 = strProteinNameRest.find(',');

		}

		outfido<<'r'<<' '<<strProteinNameRest<<endl;

		outfido<<'p'<<' '<< dProbability<<endl;


	}

}

//////////////////////////////////////////////////////////////////////////
//read from the file
void CListMassPeptideXcorr::readparsePrePareSimple_file( string strFileName )
{
	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}
	int irow=0;

	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
		string strProteinName;
		string strSpectrum;
		double dProbability;
		istringstream strStrem(buffer);
		strStrem>>strSpectrum>>strProteinName>>dProbability;
		string strPepIdentify;
		stringstream ss;
		ss<<dProbability;
		
		strPepIdentify=strSpectrum+ss.str();

		if(m_mapMassPeptidePos.find(strPepIdentify)==m_mapMassPeptidePos.end())
		{
			CMassPeptideXcorr* pMassPeptideXcorr=new CMassPeptideXcorr( strProteinName, strSpectrum, dProbability);
		
			pMassPeptideXcorr->m_vProtein.push_back(strProteinName);
			pMassPeptideXcorr->m_mapProteinPos.insert(make_pair(strProteinName,pMassPeptideXcorr->m_vProtein.size()-1));

			m_vpMassPeptideXcorr.push_back(pMassPeptideXcorr);

			m_mapMassPeptidePos.insert(make_pair(strPepIdentify,irow));

			//////////////////////////////////////////////////////////////////////////
			//insert the same peptidesequence from difference mass data
			if( m_mapvSameMassPeptides.find(strSpectrum)==m_mapvSameMassPeptides.end())
			{
				vector<int> v;
				v.push_back(irow);
				m_mapvSameMassPeptides.insert(make_pair(strSpectrum,v));
			}
			else
			{
				m_mapvSameMassPeptides[strSpectrum].push_back(irow);
			}
			irow++;



		}
		else
		{
			if(find(m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.begin(),m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.end(),strProteinName)==m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.end())
			{
				m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.push_back(strProteinName);
				m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_mapProteinPos.insert(make_pair(strProteinName,m_vpMassPeptideXcorr[m_mapMassPeptidePos[strPepIdentify]]->m_vProtein.size()-1));
			}
		}
	}
	//set average value to protein which have some shared peptide

	//////////////////////////////////////////////////////////////////////////
	//fidoFormatOutput
	ofstream outfido("fidoOutSequest.txt");
	for (int i=0;i<m_vpMassPeptideXcorr.size();i++)
	{
		outfido<<'e'<<' '<<m_vpMassPeptideXcorr[i]->m_strSpectrum<<endl;
		int iPeptideProteinSize=m_vpMassPeptideXcorr[i]->m_vProtein.size();
		for (int j = 0;j< iPeptideProteinSize;j++)
		{
			outfido<<'r'<<' '<<m_vpMassPeptideXcorr[i]->m_vProtein[j]<<endl;

		}
		outfido<<'p'<<' '<< m_vpMassPeptideXcorr[i]->m_dProbability<<endl;
	}

	SetProbPeptideToProteinOnAverage();
}

//////////////////////////////////////////////////////////////////////////
//read from the file
void CListMassPeptideXcorr::readparsePrePare_file( string strFileName )
{
	ifstream is (strFileName.c_str());
	if(!is)
	{
		throw string(("Can't open file"+strFileName).c_str());
	}
	int irow=0;

	for (string buffer;getcleanline(is,buffer) && !buffer.empty();)
	{
		string strProteinName;
		string strSpectrum;
		double dDeltaCN;
		string strPeptideSequence;
		double dMZ=0.0;//no this attribute
		double dCalculateH;
		double dActuralPeptide;
		int iSpectrumCharge;
		bool bUnique=true;
		string strFileName;
		int iScanNO;
		double dRetenTime;
		int iMissCleave;
		double dMassDiff;
		double dProbability;
		double dSp;
		bool bIPIDecoy=true;
		char cPreResidue;
		char cPostResidue;
		string strFullMod;
		int iModification;
		int iMatch;
		int iPredict;
		char cSharedUnique;
		char cIPIDecoy;

		double dXCorr=0.0;

		istringstream strStrem(buffer);//×Ö·û´®±ä³ÉÁ÷

		strStrem>>strFileName>>iScanNO>>iSpectrumCharge>>dCalculateH>>dRetenTime>>iMissCleave>>dMassDiff>>dProbability>>strProteinName>>cSharedUnique>>cIPIDecoy>>
			strPeptideSequence>>cPreResidue>>cPostResidue>>strFullMod>>iModification>>dXCorr>>dDeltaCN>>dSp>>iMatch>>iPredict;
		dActuralPeptide=dCalculateH+dMassDiff;

		if (dProbability<0.05) continue;
		
		stringstream ss;
		ss<<iScanNO;
		strSpectrum=strFileName+ss.str();//identification
		//		strStrem>>strProteinName>>strSpectrum>>dXCorr>>dDeltaCN>>strPeptideSequence>>dMZ>>dCalculateH>>dActuralPeptide>>iSpectrumCharge;
		//Peptide corresponding proteins
		if (cSharedUnique=='S')
			bUnique=false;
		if (cIPIDecoy=='D')
			bIPIDecoy=false;

		if(bUnique || m_mapMassPeptidePos.find(strSpectrum)==m_mapMassPeptidePos.end())
		{
			CMassPeptideXcorr* pMassPeptideXcorr=new CMassPeptideXcorr( strProteinName, strSpectrum, strFileName, iScanNO, dRetenTime, iMissCleave, dMassDiff, dProbability, dXCorr, dDeltaCN, dSp, strPeptideSequence, dMZ, dCalculateH, dActuralPeptide, iSpectrumCharge, bUnique, bIPIDecoy, cPreResidue, cPostResidue, strFullMod, iModification, iMatch, iPredict);
			//Protein accession numbers	Spectrum name	SEQUEST XCorr score	Peptide sequence	Observed m/z	Calculated +1H Peptide Mass (AMU)	Actual peptide mass (AMU)	Spectrum charge

			pMassPeptideXcorr->m_vProtein.push_back(strProteinName);
			pMassPeptideXcorr->m_mapProteinPos.insert(make_pair(strProteinName,pMassPeptideXcorr->m_vProtein.size()-1));

//			if (pMassPeptideXcorr->m_dProbability<0.0005) pMassPeptideXcorr->m_dProbability=0.0005;
			
			m_vpMassPeptideXcorr.push_back(pMassPeptideXcorr);

			m_mapMassPeptidePos.insert(make_pair(strSpectrum,irow));

			//////////////////////////////////////////////////////////////////////////
			//insert the same peptidesequence from difference mass data
			if( m_mapvSameMassPeptides.find(strPeptideSequence)==m_mapvSameMassPeptides.end())
			{
				vector<int> v;
				v.push_back(irow);
				m_mapvSameMassPeptides.insert(make_pair(strPeptideSequence,v));
			}
			else
			{
				m_mapvSameMassPeptides[strPeptideSequence].push_back(irow);
			}
			irow++;

			

		}
		else
		{
			m_vpMassPeptideXcorr[m_mapMassPeptidePos[strSpectrum]]->m_vProtein.push_back(strProteinName);
			m_vpMassPeptideXcorr[m_mapMassPeptidePos[strSpectrum]]->m_mapProteinPos.insert(make_pair(strProteinName,m_vpMassPeptideXcorr[m_mapMassPeptidePos[strSpectrum]]->m_vProtein.size()-1));
		}
		

	}

	//set average value to protein which have some shared peptide
	SetProbPeptideToProteinOnAverage();

	//先只考虑unique peptide
//	SetProbPeptideToProteinOnUniqueNoConsiderShared();

}

ifstream& CListMassPeptideXcorr::getcleanline( ifstream& is, string& buffer )
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

void CListMassPeptideXcorr::analyseMassPeptideXcorr()
{
	map< string,vector<int> > mapVProteinSupport;

	for (int i=0;i<m_vpMassPeptideXcorr.size();i++)
	{

		//////////////////////////////////////////////////////////////////////////
		//analyse the xcorr, the concrete parameter can be set
		string strProteinName=m_vpMassPeptideXcorr[i]->m_strProteinName;
		if(mapVProteinSupport.find(strProteinName)==mapVProteinSupport.end())
		{
			vector<int> vSupport(3,0);
			mapVProteinSupport.insert(make_pair(strProteinName,vSupport));

			
		}

		//////////////////////////////////////////////////////////////////////////
		//filter the peptides with low qualify 
		if (m_vpMassPeptideXcorr[i]->m_iSpectrumCharge==1)//one charge
		{
			if (m_vpMassPeptideXcorr[i]->m_dXcorr<2.0)
				continue;
		}else if (m_vpMassPeptideXcorr[i]->m_iSpectrumCharge==2)//tow charge
		{
			if (m_vpMassPeptideXcorr[i]->m_dXcorr<2.5)
				continue;
		}else if (m_vpMassPeptideXcorr[i]->m_iSpectrumCharge==3)//three charge
		{
			if (m_vpMassPeptideXcorr[i]->m_dXcorr<3.0)
				continue;
		}

		if (m_vpMassPeptideXcorr[i]->m_bUnique)
			mapVProteinSupport[strProteinName][1]++;//unique
		else
			mapVProteinSupport[strProteinName][2]++;//share
		mapVProteinSupport[strProteinName][0]++;//all

	}


	ofstream outf13("proteinAnalyseOut.txt");

	//////////////////////////////////////////////////////////////////////////
	//sort map and out
	vector<PROTEINPAIR> vProteinPair( mapVProteinSupport.begin(),mapVProteinSupport.end());
	sort( vProteinPair.begin(),vProteinPair.end(),cmp);
	vector< PROTEINPAIR >::iterator it;
	for(it=vProteinPair.begin();it!=vProteinPair.end();++it)
		outf13<<it->first <<'\t' << it->second[0]<< '\t' << it->second[1]<<'\t' << it->second[2]<<endl;


}
//////////////////////////////////////////////////////////////////////////
//对peptide到蛋白质进行均分
void CListMassPeptideXcorr::SetProbPeptideToProteinOnAverage()
{
	for (int i=0;i<m_vpMassPeptideXcorr.size();i++)
	{
		int iPeptideProteinSize=m_vpMassPeptideXcorr[i]->m_vProtein.size();
		for (int j = 0;j< iPeptideProteinSize;j++)
		{
			m_vpMassPeptideXcorr[i]->m_vdProbPeptideToProtein.push_back(1.0/iPeptideProteinSize);
		}
	}
}
//////////////////////////////////////////////////////////////////////////
//第一次只处理unique peptide的问题
void CListMassPeptideXcorr::SetProbPeptideToProteinOnUniqueNoConsiderShared()
{
	for (int i=0;i<m_vpMassPeptideXcorr.size();i++)
	{
		int iPeptideProteinSize=m_vpMassPeptideXcorr[i]->m_vProtein.size();
		for (int j = 0;j< iPeptideProteinSize;j++)
		{
			if (iPeptideProteinSize==1) m_vpMassPeptideXcorr[i]->m_vdProbPeptideToProtein.push_back(1.0);//对于unique peptide 第一次赋值为1
			else m_vpMassPeptideXcorr[i]->m_vdProbPeptideToProtein.push_back(0.0);//对于shared peptide 第一次赋值为0,以后在程序中调整
		}
	}
}