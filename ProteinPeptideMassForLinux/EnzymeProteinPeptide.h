#pragma once
#include <map>
#include <vector>
using namespace std;
class CEnzymeProteinPeptide
{
public:
	CEnzymeProteinPeptide(void);
	~CEnzymeProteinPeptide(void);
	map<string,string> mapEnzymeRegex;
	void SplitProteinWithEnzyme(string& strProteinName,string& strProteinSequence,string& strEnzymeType,map<string,string>& m_mapPeptideProtein);
	void SplitProteinWithEnzyme(string& strProteinSequence,string& strEnzymeType,vector<string>& vPeptides,int& ivalidPeptideCount);
};

