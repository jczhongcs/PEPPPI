#pragma once
#include <iostream>
#include <algorithm>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
/*#include <set>*/
#include <vector>
#include "Protein.h"
#include <assert.h>
// #include <assert.h>
// #include <ctime>
using namespace std;

typedef vector<int> Complex;
class CPPI
{
public:
	vector < vector<bool> > m_vvAdjcency_matrix;	//�ڽӾ���
	vector < vector<double> > m_vvWeight;		//Ȩ�ؾ���
	vector < map<int,int> > m_vmapNodeIndexForWeight;//ΪȨ���ڽӾ����node��ȡindex


// 	vector < vector<bool> > m_vvAdjcencyPG_matrix;	//�ڽӾ���
// 	vector < vector<double> > m_vvWeightPG;		//Ȩ�ؾ���


// 	vector < vector<bool> > m_vvAdjcencyGG_matrix;	//�ڽӾ���
// 	vector < vector<double> > m_vvWeightGG;		//Ȩ�ؾ���

	vector < vector<double> > m_vvgeno_express;  //����������

	vector < vector<int> > m_vvAdjcency_list;	//�ڽӱ�
	vector < vector<double> > m_vvGeno_express;  //����������
	int m_inum;					    //��������Ŀ
	map <string,int> m_mapProtein_id;	//���������Ƶ�id��ӳ���
	vector <string> m_vProtein_name;	//��������

// 	map <string,int> m_mapProteinGO_id;	//���������Ƶ�id��ӳ���
// 	vector <string> m_vProteinGO_name;	//��������
// 
// 	map <string,int> m_mapGO_id;	//���������Ƶ�id��ӳ���
// 	vector <string> m_vGO_name;	//��������

	vector <bool> m_vIs_ess;			//�Ƿ�Ϊ�ؼ�������
 	vector <Complex> m_M1;			//�ؼ����׻�õĸ�����
 	vector <Complex> m_M2;            //�ǹؼ����׹��ɵĸ�����
	vector <int> m_vtwohop;            //2����

	vector <bool> m_is_ess;			//�Ƿ�Ϊ�ؼ�������
	int m_num;					    //��������Ŀ
	double m_argR;					//��ֵ
	double m_argT;					//��ֵ
	CPPI(void);
	~CPPI(void);


	void ReadFile(string ppi_file);
	void GenerateIPI_PPI(CListProtein* plistProtein);


	double GetECC(int x,int y);

	double GetPCC(int x,int y);

	//�����Ȩ�� Ȩֵ=ECC+PCC
	void GetWeight();

	//compare complexC1 with complexC2
	double Similarity(const Complex &c1,const Complex &c2);

	//���������ɵĸ������c�ܷ���븴���Ｏ��v
	bool Test(const Complex &c, vector<Complex> &v);

	double Similarity2(const Complex &c1,const Complex &c2);

	bool Test2(const Complex &c, vector<Complex> &v);

	//��ȡ���ڹؼ������ʵĸ������
	void FindCoreEss();

	//������ڹؼ������ʵĸ������
	void ExpandEss();

	//���������
	int TwoHopDegree(int id);

	//���㸴�����ܶ�
	double Density(const Complex &c);

	//��ȡ���ڷǹؼ������ʵĸ������
	void FindCoreNoEss();

	//������ڷǹؼ������ʵĸ������
	void ExpandNoEss();
	void ReadFileWeighted(string ppi_file);
	
};

