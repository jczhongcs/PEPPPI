#pragma once
#include "calculateproteinprobability.h"
class CCalProteinProbUnique :
	public CCalculateProteinProbability
{
public:
	CCalProteinProbUnique(CListMassPeptideXcorr* pLMX,CListProtein* pLP):CCalculateProteinProbability(pLMX,pLP){};

	void calculateProteinProbability();
	
	~CCalProteinProbUnique(void);
};

