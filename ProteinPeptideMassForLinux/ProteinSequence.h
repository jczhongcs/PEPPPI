#pragma once
#include <string>
#include <vector>
#include <map>
#include "ListMassPeptideXcorr.h"
using namespace std;
class CProteinSequence
{
public:
	vector<string> m_vstrName;
	
	vector<int> m_viIsoform;
	vector< map< string, vector<string> > > m_vmapvIsoformPeptide;
	vector< string> m_vstrProteinLine;
	void parse_file(string strProteinSequenceName);

	CProteinSequence(void);
	~CProteinSequence(void);
	ifstream& getcleanline( ifstream& is, string& buffer );
	void PeptideProteinAlignment(CListMassPeptideXcorr* pListMassPeptideXcorr);
	int Alignment(  string& strSource, string& strPattern);
	int CountEnzyme(string &strPattern);
	void PeptideProteinAlignmentWithLine( CListMassPeptideXcorr* pListMassPeptideXcorr );
	void parse_fileLine( string strProteinSequenceName );

	void parse_filePeptideProteinMapp( string strfilePeptideProteinMapp );
	void parse_fileProteinIPI_Swissport_Map( string strfilePeptideProteinMapp );
	void parse_filePPIMappingProteinIPI_Swissport_Map(string strfilePPIMapp);
	void parse_PeptideProteinSequest( string strfilePPIMapp );
	void parse_ProteinIPIID_AC( string strfilePPIMapp );
};