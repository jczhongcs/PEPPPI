#pragma once
#include "calculateproteinprobability.h"
class CCalProteinProbLogisticRegression :
	public CCalculateProteinProbability
{
public:
	void calculateProteinProbability();
	CCalProteinProbLogisticRegression(CListMassPeptideXcorr* pLMX,CListProtein* pLP):CCalculateProteinProbability(pLMX,pLP){};

	~CCalProteinProbLogisticRegression(void);
};

