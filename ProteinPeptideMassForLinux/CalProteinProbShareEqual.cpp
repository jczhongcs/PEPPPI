#include "stdafx.h"
#include "CalProteinProbShareEqual.h"

#include <iostream>
#include <set>



CCalProteinProbShareEqual::~CCalProteinProbShareEqual(void)
{
}

void CCalProteinProbShareEqual::calculateProteinProbabilityOut(ofstream& out)
{
	double dProbability=1.0;
	set<string> setPeptideName;//考虑多个peptide对应一个peptide，但probability不同的情况
	for (CMassPeptideXcorr* pPeptideXcorr:m_pProtein->m_vSupportPeptide)
	{
/*
		int iPeptideSharedCountInProtein=1;
		if (pPeptideXcorr->m_bUnique==false)//shared peptide 均分
		{
			iPeptideSharedCountInProtein=pPeptideXcorr->m_vProtein.size();
		}


		dProbability=dProbability*(1-pPeptideXcorr->m_dProbability/iPeptideSharedCountInProtein);
*/

		//////////////////////////////////////////////////////////////////////////
		//test
// 		if (m_pProtein->m_strProteinName=="YOR133W")
// 			cout<<m_pProtein;
		//////////////////////////////////////////////////////////////////////////
		//filter out the shared peptide is 0, 邻居支持度为0，即shared peptide不映射该蛋白质
		if (pPeptideXcorr->m_vdProbPeptideToProtein[pPeptideXcorr->m_mapProteinPos[m_pProtein->m_strProteinName]]>0)
		{
			out<<pPeptideXcorr->m_strPeptideSequence<<'\t'<< (pPeptideXcorr->m_vdProbPeptideToProtein[pPeptideXcorr->m_mapProteinPos[m_pProtein->m_strProteinName]])<<'\t';
			dProbability=dProbability*(1-pPeptideXcorr->m_dProbability*pPeptideXcorr->m_vdProbPeptideToProtein[pPeptideXcorr->m_mapProteinPos[m_pProtein->m_strProteinName]]);
		setPeptideName.insert(pPeptideXcorr->m_strPeptideSequence);
		}
	}
	int iPeptideNameSize=1;
	if (setPeptideName.size()<m_pProtein->m_ivalidTheoryPeptideCount)
	{
		iPeptideNameSize=setPeptideName.size();
	}else
	{
		iPeptideNameSize=m_pProtein->m_ivalidTheoryPeptideCount;
	}
	m_pProtein->m_dProbability=pow((1.0-dProbability),-log((double)iPeptideNameSize/m_pProtein->m_ivalidTheoryPeptideCount));
	out<<m_pProtein->m_dProbability;
//	m_pProtein->m_dProbability=pow((1.0-dProbability),-log((double)m_pProtein->m_vSupportPeptide.size()/m_pProtein->m_ivalidTheoryPeptideCount));
//	m_pProtein->m_dProbability=(1-dProbability)*m_pProtein->m_vSupportPeptide.size()/m_pProtein->m_vTheoryPeptide.size();
}
