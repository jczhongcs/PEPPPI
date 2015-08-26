#include "stdafx.h"
#include "CalProteinProbUnique.h"




CCalProteinProbUnique::~CCalProteinProbUnique(void)
{
}

void CCalProteinProbUnique::calculateProteinProbability()
{
	double dProbability=1.0;
	for (CMassPeptideXcorr* pPeptideXcorr:m_pProtein->m_vSupportPeptide)
	{
		dProbability=dProbability*(1-pPeptideXcorr->m_dProbability);
	}
	m_pProtein->m_dProbability=pow((1-dProbability),-log(m_pProtein->m_vSupportPeptide.size()/m_pProtein->m_ivalidTheoryPeptideCount));
}
