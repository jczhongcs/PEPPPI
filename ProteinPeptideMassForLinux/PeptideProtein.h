#pragma once
#include <map>
using namespace std;
class CPeptideProtein
{
public:
	int m_iType;
	map<string,string> m_mapTheoryPeptideProtein;
	map<string,string> m_mapMassPeptideProtein;
	CPeptideProtein(void);
	~CPeptideProtein(void);
};

