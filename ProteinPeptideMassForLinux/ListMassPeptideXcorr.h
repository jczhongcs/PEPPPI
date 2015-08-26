#pragma once
#include <string>
#include <vector>
#include <map>
//#include "Protein.h"
using namespace std;
class CMassPeptideXcorr
{
public:
	string m_strProteinName;//蛋白质名称
	string m_strSpectrum;//
	
	string m_strFileName;
	int	   m_iScanNo;//扫描号
	double m_dRetentionTime;//保留时间
	int    m_iMissCleave;//错过切开数
	double m_dMassDiff;//与理论谱的质量差
	double m_dProbability;//存在概率
	
	double m_dXcorr;//sequest xcorr
	double m_dDeltaCN;//sequest deltaCN
	double m_dSp;//sequest SP

	string m_strPeptideSequence;//肽序列
	double m_dMZ;//质荷比
	double m_dCalculateH;//
	double m_dActuralPeptide;
	int m_iSpectrumCharge;

	bool m_bUnique;//是否唯一
	bool m_bIPI_Decoy;//是不是IPI—Decoy
	char m_cPreResidue;//preResidue
	char m_cPostResidue;//PostResidue

	string m_strFullModofication;//全修饰
	int m_iModificationCount;//修饰个数

	int m_iMatchingIon;
	int m_iPredictedIon;

	vector<string> m_vProtein;//映射的蛋白质，若为unique则size为1
	vector<double> m_vdProbPeptideToProtein;//映射到蛋白质的概率，若为unique则size为1，且值为1.0，否则size>1,而累加值为1.0;
	map<string ,int> m_mapProteinPos;//protein的位置
	void CalculateProbPeptideToProtein();
/*
	string m_strName;
	string m_strPeptide;*/
	

	CMassPeptideXcorr(string strProteinName,string strSpectrum,double dXcorr,double dDeltaCN,string strPeptideSequence,double dMZ,double dCalculateH,double dActuralPeptide,int iSpectrumCharge,bool bUnique): m_strProteinName(strProteinName), m_strSpectrum(strSpectrum),m_dXcorr(dXcorr),m_dDeltaCN(dDeltaCN), m_strPeptideSequence(strPeptideSequence), m_dMZ(dMZ), m_dCalculateH(dCalculateH), m_dActuralPeptide(dActuralPeptide), m_iSpectrumCharge(iSpectrumCharge), m_bUnique(bUnique){}

	CMassPeptideXcorr(string strProteinName,string strSpectrum,string strFileName,int iScanNO,double dRetenTime,int iMissCleave,double dMassDiff,double dProbability,double dXcorr,double dDeltaCN,double dSp,string strPeptideSequence,double dMZ,double dCalculateH,double dActuralPeptide,int iSpectrumCharge,bool bUnique,bool bIPIDecoy,char cPreResidue,char cPostResidue,string strFullMod,int iModificationCount,int iMatch,int iPredict): m_strProteinName(strProteinName), m_strSpectrum(strSpectrum),m_strFileName(strFileName),m_iScanNo(iScanNO),m_dRetentionTime(dRetenTime),m_iMissCleave(iMissCleave),m_dMassDiff(dMassDiff), m_dProbability(dProbability),m_dXcorr(dXcorr),m_dDeltaCN(dDeltaCN),m_dSp(dSp), m_strPeptideSequence(strPeptideSequence), m_dMZ(dMZ), m_dCalculateH(dCalculateH), m_dActuralPeptide(dActuralPeptide), m_iSpectrumCharge(iSpectrumCharge), m_bUnique(bUnique),m_bIPI_Decoy(bIPIDecoy),m_cPreResidue(cPreResidue),m_cPostResidue(cPostResidue),
		m_strFullModofication(strFullMod),m_iModificationCount(iModificationCount),m_iMatchingIon(iMatch),m_iPredictedIon(iPredict){}
	CMassPeptideXcorr(string strProteinName,string strSpectrum,double dProb):m_strProteinName(strProteinName), m_strSpectrum(strSpectrum),m_dProbability(dProb),m_strPeptideSequence(strSpectrum){}
};
class CListMassPeptideXcorr
{
public:
	CListMassPeptideXcorr(void);
	~CListMassPeptideXcorr(void);
	vector< CMassPeptideXcorr* > m_vpMassPeptideXcorr;//肽列表
	map< string, int > m_mapMassPeptidePos;//肽列表位置

	map< string, vector<int>> m_mapvSameMassPeptides;//save the same peptidesequence from difference mass data


	int m_iType;
	void readparse_file(string strFileName,int iType=1/*sequest*/);
	ifstream& getcleanline( ifstream& is, string& buffer );
	void analyseMassPeptideXcorr();
	void readparsePrePare_file( string strFileName );

	void SetProbPeptideToProteinOnAverage();

	void readparsePrePareSimple_file( string strFileName );
	void readparsePrePareSimple_file( string strFileName,char charSplitChar );
	void readparseFidoe_file( string strFileName );
	void readparseFidoe_file( string strFileName,char charSplitChar );
	void convertToFido( string strFileName );
	void SetProbPeptideToProteinOnUniqueNoConsiderShared();
};

