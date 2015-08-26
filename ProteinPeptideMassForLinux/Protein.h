#pragma once
#include <string>
#include <map>
#include "PeptideProtein.h"
#include <vector>
#include <set>
#include "ListMassPeptideXcorr.h"
#include "EnzymeProteinPeptide.h"
using namespace std;
class CProtein
{
public:

	string m_strProteinName;//蛋白质名称
	vector<string> m_strSwissName;//SwissName s
	string m_strIPIName;//IPI
	vector<string> m_strTREMBL;//trembl
	vector<string> m_strENSEMBL;//ensembl
	vector<string> m_strREFSEQ;//refseq

	string m_strSequence;//sequence
	string m_strHead;//head
	string m_strGeneName;//GeneName..yeast
	CProtein(string& strHead,string& strSequence);
	CProtein(string& strHead,string& strSequence,bool bYeast);
	CProtein(string& strHead,string& strSequence,bool bHuman,bool bWithDecoy);
	CProtein(string& strHead,string& strSequence,bool bHuman,bool bWithDecoy,bool ensembl);

	void ExtractNames( string &strName, vector<string>& vstr );

	~CProtein(void);
//	vector< map<string,int> > m_vSupportPeptide;//peptide isShare
	vector< CMassPeptideXcorr* > m_vSupportPeptide;//支持的肽

	double m_dSupportValueFromOtherProtein;//来源于邻居的支持度

	vector<string> m_vTheoryPeptide;//理论肽
	int m_ivalidTheoryPeptideCount;//有效理论肽个数
	double m_dProbability;//蛋白质概率
	//读取三个输入文件
	CEnzymeProteinPeptide objEnzym;//剪切enzyme
	bool m_bDecoy;

	void parse_Sequence(CPeptideProtein& peptideProtein);
	bool findProtein(string str);//根据str(m_strSwissName|m_strTREMBL|m_strENSEMBL|m_strREFSEQ)查找蛋白质
	void outFileProteinName(ofstream& out);//输出蛋白质名称集
	void outFilePeptideProbability( ofstream& out,double myValue=0.99 ,bool bShowDetail=false);//根据double myValue输出自信peptide，并输出全部支持Peptide


	//////////////////////////////////////////////////////////////////////////
	//calculate the score with shared peptide, which view as unique peptide.
	void calculateProbabilitySharedOne();
	//////////////////////////////////////////////////////////////////////////
	//calculate the score with shared peptide, which divide equally into protein  
	void calculateProbabilitySharedDivide();
	//////////////////////////////////////////////////////////////////////////
	//calculate the score with shared peptide, which is calculated the probability with logistic regression 
	void calculateProbabilitySharedLogisticRegression();
	//////////////////////////////////////////////////////////////////////////
	//calculate the score with shared peptide, which is calculated the probability with wu's methods 
//	void calculateProbabilitySharedWeighted(CListMassPeptideXcorr* pListMassPeptideXcorr,CListProtein* pListProtein);

	double getSameMassSPValue( CListMassPeptideXcorr* pListMassPeptideXcorr, CMassPeptideXcorr* pPeptideXcorr );
	bool findProteinSequence( string str );


};


class CListProtein
{
public:
	CListProtein(void){};
	~CListProtein(void);
	map< string ,int > mapListPos;//蛋白质在vProteinList的位置
	vector<CProtein* > vProteinList;//蛋白质列表
	void readparse_file(string strFileName);//从文件中读取蛋白质列表
	void OutFileProteinListHeadInfo(string strFileName);//输入蛋白质列表中的蛋白质名称
	void MatchProteinHeadList(string str,vector<string>& vIPIName);//查找蛋白质列表中蛋白质名称为str的蛋白质
	ifstream& getcleanline( ifstream& is, string& buffer );
	void mergePPI( string str );//合并PPI
	vector<CMassPeptideXcorr*> getSupportPeptideFromProteinStr( string strProtein )//获取指定蛋白质的支持肽
	{
		return vProteinList[mapListPos[strProtein]]->m_vSupportPeptide;
	}
	void readparse_Yeastfile( string strFileName );
	bool isUniquePeptide(CMassPeptideXcorr* pPeptide);
	bool isPeptideInProtein(CMassPeptideXcorr* pPeptide,ofstream& out);
	void outProteinAndSequence( CProtein* objProtein ,int iCount=100);

	void readGoldenSet(string strFileName);

	vector<string> vYeastGoldenSet;//蛋白质参考集

	int isInGoldenset(string proteinName);

	int iGetProteinPos(string proteinName);
	bool OutNoRepeatPeptideInProtein( CMassPeptideXcorr* pPeptide,ofstream& out );
	void readparse_Esembelfile( string strFileName );
};
