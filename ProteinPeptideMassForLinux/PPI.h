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
	vector < vector<bool> > m_vvAdjcency_matrix;	//邻接矩阵
	vector < vector<double> > m_vvWeight;		//权重矩阵
	vector < map<int,int> > m_vmapNodeIndexForWeight;//为权重邻接矩阵的node获取index


// 	vector < vector<bool> > m_vvAdjcencyPG_matrix;	//邻接矩阵
// 	vector < vector<double> > m_vvWeightPG;		//权重矩阵


// 	vector < vector<bool> > m_vvAdjcencyGG_matrix;	//邻接矩阵
// 	vector < vector<double> > m_vvWeightGG;		//权重矩阵

	vector < vector<double> > m_vvgeno_express;  //基因表达数据

	vector < vector<int> > m_vvAdjcency_list;	//邻接表
	vector < vector<double> > m_vvGeno_express;  //基因表达数据
	int m_inum;					    //蛋白质数目
	map <string,int> m_mapProtein_id;	//蛋白质名称到id的映射表
	vector <string> m_vProtein_name;	//蛋白质名

// 	map <string,int> m_mapProteinGO_id;	//蛋白质名称到id的映射表
// 	vector <string> m_vProteinGO_name;	//蛋白质名
// 
// 	map <string,int> m_mapGO_id;	//蛋白质名称到id的映射表
// 	vector <string> m_vGO_name;	//蛋白质名

	vector <bool> m_vIs_ess;			//是否为关键蛋白质
 	vector <Complex> m_M1;			//关键蛋白获得的复合物
 	vector <Complex> m_M2;            //非关键蛋白构成的复合物
	vector <int> m_vtwohop;            //2极度

	vector <bool> m_is_ess;			//是否为关键蛋白质
	int m_num;					    //蛋白质数目
	double m_argR;					//阈值
	double m_argT;					//阈值
	CPPI(void);
	~CPPI(void);


	void ReadFile(string ppi_file);
	void GenerateIPI_PPI(CListProtein* plistProtein);


	double GetECC(int x,int y);

	double GetPCC(int x,int y);

	//计算边权重 权值=ECC+PCC
	void GetWeight();

	//compare complexC1 with complexC2
	double Similarity(const Complex &c1,const Complex &c2);

	//测试新生成的复合物核c能否加入复合物集合v
	bool Test(const Complex &c, vector<Complex> &v);

	double Similarity2(const Complex &c1,const Complex &c2);

	bool Test2(const Complex &c, vector<Complex> &v);

	//获取基于关键蛋白质的复合物核
	void FindCoreEss();

	//扩充基于关键蛋白质的复合物核
	void ExpandEss();

	//计算二级度
	int TwoHopDegree(int id);

	//计算复合物密度
	double Density(const Complex &c);

	//获取基于非关键蛋白质的复合物核
	void FindCoreNoEss();

	//扩充基于非关键蛋白质的复合物核
	void ExpandNoEss();
	void ReadFileWeighted(string ppi_file);
	
};

