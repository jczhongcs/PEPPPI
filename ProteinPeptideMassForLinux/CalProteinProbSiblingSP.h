#pragma once
#include "calculateproteinprobability.h"
class CCalProteinProbSiblingSP :
	public CCalculateProteinProbability
{
public:
	CCalProteinProbSiblingSP(CListMassPeptideXcorr* pLMX,CListProtein* pLP):CCalculateProteinProbability(pLMX,pLP){};
	void calculateProteinProbability();
	~CCalProteinProbSiblingSP(void);
};

