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

	string m_strProteinName;//����������
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
	vector< CMassPeptideXcorr* > m_vSupportPeptide;//֧�ֵ���

	double m_dSupportValueFromOtherProtein;//��Դ���ھӵ�֧�ֶ�

	vector<string> m_vTheoryPeptide;//������
	int m_ivalidTheoryPeptideCount;//��Ч�����ĸ���
	double m_dProbability;//�����ʸ���
	//��ȡ���������ļ�
	CEnzymeProteinPeptide objEnzym;//����enzyme
	bool m_bDecoy;

	void parse_Sequence(CPeptideProtein& peptideProtein);
	bool findProtein(string str);//����str(m_strSwissName|m_strTREMBL|m_strENSEMBL|m_strREFSEQ)���ҵ�����
	void outFileProteinName(ofstream& out);//������������Ƽ�
	void outFilePeptideProbability( ofstream& out,double myValue=0.99 ,bool bShowDetail=false);//����double myValue�������peptide�������ȫ��֧��Peptide


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
	map< string ,int > mapListPos;//��������vProteinList��λ��
	vector<CProtein* > vProteinList;//�������б�
	void readparse_file(string strFileName);//���ļ��ж�ȡ�������б�
	void OutFileProteinListHeadInfo(string strFileName);//���뵰�����б��еĵ���������
	void MatchProteinHeadList(string str,vector<string>& vIPIName);//���ҵ������б��е���������Ϊstr�ĵ�����
	ifstream& getcleanline( ifstream& is, string& buffer );
	void mergePPI( string str );//�ϲ�PPI
	vector<CMassPeptideXcorr*> getSupportPeptideFromProteinStr( string strProtein )//��ȡָ�������ʵ�֧����
	{
		return vProteinList[mapListPos[strProtein]]->m_vSupportPeptide;
	}
	void readparse_Yeastfile( string strFileName );
	bool isUniquePeptide(CMassPeptideXcorr* pPeptide);
	bool isPeptideInProtein(CMassPeptideXcorr* pPeptide,ofstream& out);
	void outProteinAndSequence( CProtein* objProtein ,int iCount=100);

	void readGoldenSet(string strFileName);

	vector<string> vYeastGoldenSet;//�����ʲο���

	int isInGoldenset(string proteinName);

	int iGetProteinPos(string proteinName);
	bool OutNoRepeatPeptideInProtein( CMassPeptideXcorr* pPeptide,ofstream& out );
	void readparse_Esembelfile( string strFileName );
};
