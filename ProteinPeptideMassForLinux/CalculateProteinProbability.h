#pragma once
#include "ListMassPeptideXcorr.h"
#include "Protein.h"
class CCalculateProteinProbability
{
public:
	virtual void calculateProteinProbability()=0;
	CListMassPeptideXcorr* m_pListMassPeptideXcorr;
	CListProtein* m_pListProtein;
	CProtein* m_pProtein;
	CCalculateProteinProbability(CListMassPeptideXcorr* pLMX,CListProtein* pLP):m_pListMassPeptideXcorr(pLMX),m_pListProtein(pLP){};
	~CCalculateProteinProbability(void);


	double getSameMassSPValue( CMassPeptideXcorr* pPeptideXcorr );
	virtual void calculateProteinProbabilityOut( ofstream& outTimes )=0;
	int m_iType;

};

