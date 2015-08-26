#pragma once
#include <fstream>
#include "calculateproteinprobability.h"
class CCalProteinProbShareEqual :
	public CCalculateProteinProbability
{
public:
	void calculateProteinProbabilityOut(ofstream& out);
	void calculateProteinProbability(){};
	CCalProteinProbShareEqual(CListMassPeptideXcorr* pLMX,CListProtein* pLP):CCalculateProteinProbability(pLMX,pLP){};
	
	~CCalProteinProbShareEqual(void);
};

