// ProteinPeptideMassForLinux.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "ProteinSequence.h"
#include "PPI.h"
#include "Protein.h"
#include "CalProteinProbShareEqual.h"
#include "Config.h"
// #include <Eigen/Eigen>
// #include <Eigen/SparseQR>



ifstream& getcleanline( ifstream& is, string& buffer )
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




void outGLNeiborProbCount( CPPI &objPPI, string strIPIName, ofstream& out, CListProtein &objListProtein,double myValue )
{
	int iCount=0;
	for(int j:objPPI.m_vvAdjcency_list[objPPI.m_mapProtein_id[strIPIName]])
	{
		
		if (objListProtein.vProteinList[objListProtein.mapListPos[objPPI.m_vProtein_name[j]]]->m_dProbability>=myValue)
		{
			iCount++;
		}
	}
	out<<iCount<<'\t';
}

int _tmain(int argc, _TCHAR* argv[])
{



//////////////////////////////////////////////////////////////////////////

	CProteinSequence objProteinSequence;


	CPPI objPPI;

	objPPI.ReadFile(Config::strPPIFileName);//yeast data04

	CListProtein objListProtein;


	objListProtein.readparse_Yeastfile(Config::strSequenceFileName);//yeast
	objListProtein.readGoldenSet(Config::strGoldensetFileName);//读入参考集 for data02




	CListMassPeptideXcorr objListMassPeptideXcorr;


	
	objListMassPeptideXcorr.readparsePrePareSimple_file(Config::strPeptideProphetResult,Config::cSplitChar);



	//////////////////////////////////////////////////////////////////////////
	//a protein include some supported peptides
	for (int i=0;i<objListMassPeptideXcorr.m_vpMassPeptideXcorr.size();i++)
	{
		for (int j = 0;j< objListMassPeptideXcorr.m_vpMassPeptideXcorr[i]->m_vProtein.size();j++)
		{
			if (objListProtein.mapListPos.find( objListMassPeptideXcorr.m_vpMassPeptideXcorr[i]->m_vProtein[j])!=objListProtein.mapListPos.end())
			{
		
				objListProtein.vProteinList[objListProtein.mapListPos[ objListMassPeptideXcorr.m_vpMassPeptideXcorr[i]->m_vProtein[j]]]->m_vSupportPeptide.push_back(objListMassPeptideXcorr.m_vpMassPeptideXcorr[i]);//add peptide to protein
			}
		}
		
	}

	//////////////////////////////////////////////////////////////////////////
	

	ofstream outProteinTheoryPeptideSize("outProteinTheoryPeptideSize");
	for (int i =0;i<objListProtein.vProteinList.size();i++)
	{
		outProteinTheoryPeptideSize<<objListProtein.vProteinList[i]->m_strProteinName<<'\t'<<objListProtein.vProteinList[i]->m_vTheoryPeptide.size()<<'\t'<<objListProtein.vProteinList[i]->m_ivalidTheoryPeptideCount<<endl;
	}



	for (int iTimes=0;iTimes<500;iTimes++)
	{
		

	//////////////////////////////////////////////////////////////////////////
	//calculate the initialization of probability of protein via shared peptide, which assigned equally to every protein

			ofstream outTimes("outTimes");
			outTimes<<iTimes<<endl;
	CCalculateProteinProbability* pCalculateProteinProbabilityTools=new CCalProteinProbShareEqual(&objListMassPeptideXcorr,&objListProtein);
	for (int i =0;i<objListProtein.vProteinList.size();i++)
	{
		pCalculateProteinProbabilityTools->m_pProtein=objListProtein.vProteinList[i];
		outTimes<<pCalculateProteinProbabilityTools->m_pProtein->m_strProteinName<<"\t";
		pCalculateProteinProbabilityTools->calculateProteinProbabilityOut(outTimes);
		outTimes<<endl;
	}





	//////////////////////////////////////////////////////////////////////////
	//calculate supportScore
	for (int i =0;i<objListProtein.vProteinList.size();i++)
	{
		double dNeiborProbability=0.0;
//		double dDifferWithNeiborProbability=0.0;

		if (objPPI.m_mapProtein_id.find(objListProtein.vProteinList[i]->m_strIPIName)!=objPPI.m_mapProtein_id.end())//mapping to PPI
		{
			int iNeiborSupportSize=0;
			int iNeiborHightQualitySize=0;
			int iNeiborPoorQualitySize=0;
			int iProteinMapid=objPPI.m_mapProtein_id[objListProtein.vProteinList[i]->m_strIPIName];
			for(int j:objPPI.m_vvAdjcency_list[iProteinMapid])
			{
				
				if (objListProtein.vProteinList[objListProtein.mapListPos[objPPI.m_vProtein_name[j]]]->m_dProbability<0.1)
				{
					iNeiborPoorQualitySize++;

				}else if (objListProtein.vProteinList[objListProtein.mapListPos[objPPI.m_vProtein_name[j]]]->m_dProbability>0.9)
				{
					iNeiborHightQualitySize++;
				}
				
				
				if (objListProtein.vProteinList[objListProtein.mapListPos[objPPI.m_vProtein_name[j]]]->m_dProbability*objPPI.m_vvWeight[iProteinMapid][objPPI.m_vmapNodeIndexForWeight[iProteinMapid][j]]>0.0)
				{
				
					
					dNeiborProbability+=objListProtein.vProteinList[objListProtein.mapListPos[objPPI.m_vProtein_name[j]]]->m_dProbability*objPPI.m_vvWeight[iProteinMapid][objPPI.m_vmapNodeIndexForWeight[iProteinMapid][j]];	
					
					iNeiborSupportSize++;
				}
			}
			


			if(iNeiborSupportSize)
				dNeiborProbability/=iNeiborSupportSize;
			else
				dNeiborProbability=0.0;

		}else//no mapping to PPI
		{
			dNeiborProbability=0.0;
		}
		objListProtein.vProteinList[i]->m_dSupportValueFromOtherProtein=dNeiborProbability;
//		objListProtein.vProteinList[i]->m_dProbability=(objListProtein.vProteinList[i]->m_dProbability+dDifferWithNeiborProbability);
	}

	//////////////////////////////////////////////////////////////////////////
	//adjust shared Peptide
	for (CMassPeptideXcorr* pMassPepXcorr:objListMassPeptideXcorr.m_vpMassPeptideXcorr)
	{
		if (pMassPepXcorr->m_vProtein.size()>1)
		{
		
		
			double dProbSupportAndProtein=0.0;
			vector<double> vProbPeptidetoProtein(pMassPepXcorr->m_vProtein.size(),0.0);

			for (int j = 0;j< pMassPepXcorr->m_vProtein.size();j++)//peptide 
			{
				vProbPeptidetoProtein[j]=objListProtein.vProteinList[objListProtein.mapListPos[pMassPepXcorr->m_vProtein[j]]]->m_dSupportValueFromOtherProtein * objListProtein.vProteinList[objListProtein.mapListPos[pMassPepXcorr->m_vProtein[j]]]->m_dProbability;
				dProbSupportAndProtein+=vProbPeptidetoProtein[j];
			}
			if (dProbSupportAndProtein>0.0)
				for (int j = 0;j< pMassPepXcorr->m_vProtein.size();j++)
				{
					
					pMassPepXcorr->m_vdProbPeptideToProtein[j]=min(pMassPepXcorr->m_vdProbPeptideToProtein[j]+vProbPeptidetoProtein[j]/dProbSupportAndProtein,1.0);
				}			
		

		}
	}

	}//end times
	ofstream out(Config::strOutResultFileName);
	int irrevanceCount=0;
	int iUnique=0;
	for (int i=0;i<objListProtein.vProteinList.size();i++)
	{

		//////////////////////////////////////////////////////////////////////////
		//get # shared and unique peptide from a protein 
		int iNoUniquePep,iNoSharedPep;
		iNoUniquePep=iNoSharedPep=0;
		for (CMassPeptideXcorr* pPeptideXcorr:objListProtein.vProteinList[i]->m_vSupportPeptide)
		{

		}
		//////////////////////////////////////////////////////////////////////////
		//mapping to PPI
		string strIPIName=objListProtein.vProteinList[i]->m_strIPIName;
		if (objPPI.m_mapProtein_id.find(strIPIName)!=objPPI.m_mapProtein_id.end())
		{


			if (objListProtein.vProteinList[i]->m_vSupportPeptide.size()>0)
			{
		
				out<<strIPIName<<"\t"<<objListProtein.vProteinList[i]->m_dProbability<<'\t'<<objListProtein.isInGoldenset(objListProtein.vProteinList[i]->m_strIPIName)<<'\t'<<objListProtein.vProteinList[i]->m_vSupportPeptide.size()<<"\tNei:"<<objPPI.m_vvAdjcency_list[objPPI.m_mapProtein_id[strIPIName]].size()<<'\t';

					for (int iPercent=10;iPercent>0;iPercent--)
				{
						outGLNeiborProbCount(objPPI, strIPIName, out, objListProtein,(double)iPercent/10);

				}
						out<<endl;
			}
			
		}else//no mapping to PPI
		{
			if (objListProtein.vProteinList[i]->m_vSupportPeptide.size()>0)
			{
			out<<strIPIName<<"\t"<<objListProtein.vProteinList[i]->m_dProbability<<'\t'<<objListProtein.isInGoldenset(objListProtein.vProteinList[i]->m_strIPIName)<<'\t'<<objListProtein.vProteinList[i]->m_vSupportPeptide.size()<<"\tnoMapping\t";
					out<<endl;
			}
	
		}

	}
	cout<<irrevanceCount<<endl;
	cout<<iUnique<<endl;



	return 0;
}

