#include "stdafx.h"
#include "CalculateProteinProbability.h"




CCalculateProteinProbability::~CCalculateProteinProbability(void)
{
}

double CCalculateProteinProbability::getSameMassSPValue( CMassPeptideXcorr* pPeptideXcorr )
{
	double dPeptideSP=0.0;
	for (int i=0;i<m_pListMassPeptideXcorr->m_mapvSameMassPeptides[pPeptideXcorr->m_strPeptideSequence].size();i++)
	{
		dPeptideSP+=m_pListMassPeptideXcorr->m_vpMassPeptideXcorr[m_pListMassPeptideXcorr->m_mapvSameMassPeptides[pPeptideXcorr->m_strPeptideSequence][i]]->m_dSp;
	}	return dPeptideSP;
}

