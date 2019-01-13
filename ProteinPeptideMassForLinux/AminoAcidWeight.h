#pragma once

#include <map>
#include <string>
//#include <regex>

using namespace std;
class CAminoAcidWeight
{
public:
	CAminoAcidWeight();
	~CAminoAcidWeight();

	map<string,double> m_mapMonoisotopicWeight;
	map<string,double> m_mapAverageWeight;

	const static double m_rc[26];

	const static double m_rcnt[26];

	const static double m_nt[3];

	double getMonoisotopicWeight(string& strPeptideSequence);
	double getAverageWeight(string& strPeptideSequence);

	double hydrophobicity(string& bytes, int start, int n);

private:
	const static double m_minMass;
	const static double m_maxMass;
};

