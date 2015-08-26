#include "stdafx.h"
#include "PPI.h"



CPPI::CPPI(void)
{
}


CPPI::~CPPI(void)
{
}




void CPPI::ReadFile(string ppi_file)
{
	//读入ppi网络,首先对蛋白质进行标号
	ifstream ppi(ppi_file.c_str());
	string s1,s2;
	int index=0;
	while(ppi>>s1)
	{
		if(m_mapProtein_id.find(s1)==m_mapProtein_id.end())
		{
			m_mapProtein_id.insert(make_pair(s1,index++));
			m_vProtein_name.push_back(s1);
			map<int,int> mapNodeIndex;
//			m_vvWeight.push_back(vector<double>());
			m_vmapNodeIndexForWeight.push_back(mapNodeIndex);
		}
	}
	m_inum=index;

	m_vIs_ess.resize(m_inum);
//	m_vvAdjcency_matrix.resize(m_inum);
	m_vvAdjcency_list.resize(m_inum);
	m_vvWeight.resize(m_inum);
	
	m_vvGeno_express.resize(m_inum);
/*

	for(int i=0;i<m_inum;++i)
	{
//		m_vvAdjcency_matrix[i].resize(m_inum);
		m_vvWeight[i].resize(m_inum);
	}*/
	ppi.close();
	ppi.open(ppi_file.c_str(),ifstream::in);
	while(ppi>>s1>>s2)
	{
		int id1=m_mapProtein_id[s1];
		int id2=m_mapProtein_id[s2];
//		m_vvAdjcency_matrix[id1][id2]=1;
//		m_vvAdjcency_matrix[id2][id1]=1;
		m_vvWeight[id1].push_back(1.0);
		m_vmapNodeIndexForWeight[id1][id2]=m_vmapNodeIndexForWeight[id1].size();
//		m_vvWeight[id1][id2]=1.0;
		m_vvWeight[id2].push_back(1.0);
		m_vmapNodeIndexForWeight[id2][id1]=m_vmapNodeIndexForWeight[id2].size();
//		m_vvWeight[id2][id1]=1.0;
		m_vvAdjcency_list[id1].push_back(id2);
		m_vvAdjcency_list[id2].push_back(id1);
	}
	ppi.close();

/*

	//读取关键蛋白质文件
	ifstream ess(ess_file.c_str());
	while(ess>>s1)
	{
		if(protein_id.find(s1)!=protein_id.end())
			is_ess[protein_id[s1]]=1;
	}
	ess.close();*/

}



void CPPI::ReadFileWeighted(string ppi_file)
{
	//读入ppi网络,首先对蛋白质进行标号
	ifstream ppi(ppi_file.c_str());
	string s1,s2;
	int index=0;
	int iPos=0;
	while(ppi>>s1)
	{
		iPos++;
		if (iPos==3)
		{
			iPos=0;
		}else
		{
			if(m_mapProtein_id.find(s1)==m_mapProtein_id.end())
			{
				m_mapProtein_id.insert(make_pair(s1,index++));
				m_vProtein_name.push_back(s1);
			}
		}
	}
	m_inum=index;
	m_vIs_ess.resize(m_inum);
//	m_vvAdjcency_matrix.resize(m_inum);
	m_vvAdjcency_list.resize(m_inum);
	m_vvWeight.resize(m_inum);
	m_vvGeno_express.resize(m_inum);

	for(int i=0;i<m_inum;++i)
	{
//		m_vvAdjcency_matrix[i].resize(m_inum);
		m_vvWeight[i].resize(m_inum);
	}
	ppi.close();
	ppi.open(ppi_file.c_str(),ifstream::in);
	while(ppi>>s1>>s2)
	{
		double dWeighted=0.0;
		ppi>>dWeighted;
		int id1=m_mapProtein_id[s1];
		int id2=m_mapProtein_id[s2];
//		m_vvAdjcency_matrix[id1][id2]=1;
//		m_vvAdjcency_matrix[id2][id1]=1;

		m_vvWeight[id1].push_back(dWeighted);
		m_vmapNodeIndexForWeight[id1][id2]=m_vmapNodeIndexForWeight[id1].size();
		//		m_vvWeight[id1][id2]=1.0;
		m_vvWeight[id2].push_back(dWeighted);
		m_vmapNodeIndexForWeight[id2][id1]=m_vmapNodeIndexForWeight[id2].size();
// 		m_vvWeight[id1][id2]=dWeighted;
// 		m_vvWeight[id2][id1]=dWeighted;
		m_vvAdjcency_list[id1].push_back(id2);
		m_vvAdjcency_list[id2].push_back(id1);
	}
	ppi.close();

/*

	//读取关键蛋白质文件
	ifstream ess(ess_file.c_str());
	while(ess>>s1)
	{
		if(protein_id.find(s1)!=protein_id.end())
			is_ess[protein_id[s1]]=1;
	}
	ess.close();*/

}

void CPPI::GenerateIPI_PPI( CListProtein* plistProtein )
{
	map<string,int> mapIPIProtein_id;

	vector< vector<int> > vvIPI_PPI;
	vector<string> vIPIName;

	vector< string> vIPINameP1;
	vector< string> vIPINameP2;
	int index=0;
	for (int i =0 ;i<m_vvAdjcency_list.size();i++)
	{
		vIPINameP1.clear();
		plistProtein->MatchProteinHeadList(m_vProtein_name[i],vIPINameP1);
		for (int j=0;j<m_vvAdjcency_list[i].size();j++)
		{
			vIPINameP2.clear();
			plistProtein->MatchProteinHeadList(m_vProtein_name[m_vvAdjcency_list[i][j]],vIPINameP2);
			
			for (int i1=0;i1<vIPINameP1.size();i1++)
			{
				if(mapIPIProtein_id.find(vIPINameP1[i1])==mapIPIProtein_id.end())
				{
						mapIPIProtein_id.insert(make_pair(vIPINameP1[i1],index++));
						vIPIName.push_back(vIPINameP1[i1]);
					}
				
			}
			for (int j1=0;j1<vIPINameP2.size();j1++)
			{
				if(mapIPIProtein_id.find(vIPINameP2[j1])==mapIPIProtein_id.end())
				{
					mapIPIProtein_id.insert(make_pair(vIPINameP2[j1],index++));
					vIPIName.push_back(vIPINameP2[j1]);
				}
			}
			vvIPI_PPI.resize(index);
			for (int i1=0;i1<vIPINameP1.size();i1++)
			{
				for (int j1=0;j1<vIPINameP2.size();j1++)
				{

						int id1=mapIPIProtein_id[vIPINameP1[i1]];
						int id2=mapIPIProtein_id[vIPINameP2[j1]];
//						m_vvAdjcency_matrix[id1][id2]=1;
//						m_vvAdjcency_matrix[id2][id1]=1;
						if (find(vvIPI_PPI[id1].begin(),vvIPI_PPI[id1].end(),id2)==vvIPI_PPI[id1].end())
							vvIPI_PPI[id1].push_back(id2);
						if (find(vvIPI_PPI[id2].begin(),vvIPI_PPI[id2].end(),id1)==vvIPI_PPI[id2].end())
							vvIPI_PPI[id2].push_back(id1);
					
				}
			}
		}

	}


	ofstream out13("IPIList");
	for (int i=0;i<vIPIName.size();i++)
	{
		out13<<i<<'\t'<<vIPIName[i]<<endl;
	}
	ofstream out14("mapIPIProtein_id");
	map<string,int>::iterator   it=mapIPIProtein_id.begin();   
	for(;it!=mapIPIProtein_id.end();++it)   
		out14<<it->first<<'\t'<<it->second<<endl;   

	ofstream out15("IPIPPI");
	for (int i1=0;i1<vvIPI_PPI.size();i1++)
	{
		for (int j1=0;j1<vvIPI_PPI[i1].size();j1++)
		{

			out15<<i1<<'\t'<<vvIPI_PPI[i1][j1]<<endl;

		}
	}
	
}

double CPPI::GetECC( int x,int y )
{
	int kx=m_vvAdjcency_list[x].size();
	int ky=m_vvAdjcency_list[y].size();
	if(kx<=1||ky<=1)
		return 0;
	int Z3=0;
	for(int n1:m_vvAdjcency_list[x])
		for(int n2:m_vvAdjcency_list[y])
			if(n1==n2) ++Z3;
	return (double)Z3/min(kx-1,ky-1);
}

double CPPI::GetPCC( int x,int y )
{
	if(m_vvgeno_express[x].empty()||m_vvgeno_express[y].empty())
		return -1;
	assert(m_vvgeno_express[x].size()==m_vvgeno_express[y].size());
	int N=m_vvgeno_express[x].size();
	double xa=0,ya=0;
	for(double d:m_vvgeno_express[x])
		xa+=d;
	for(double d:m_vvgeno_express[y])
		ya+=d;
	xa/=N;
	ya/=N;
	double t1=0,t2=0,t3=0;
	for(int i=0;i<N;++i)
	{
		double xi=m_vvgeno_express[x][i];
		double yi=m_vvgeno_express[y][i];
		t1+=(xi-xa)*(yi-ya);
		t2+=(xi-xa)*(xi-xa);
		t3+=(yi-ya)*(yi-ya);
	}
	return t1/sqrt(t2*t3);
}

void CPPI::GetWeight()
{
	ofstream os("debug_weight.txt");

	os<<"A\t B\t ECC\t PCC\t Weight"<<endl;
	for(int i=0;i<m_num;++i)
		for(int j=i+1;j<m_num;++j)
		{

			if(m_vvAdjcency_matrix[i][j])
			{
				m_vvWeight[j][i]=m_vvWeight[i][j]=GetECC(i,j)+GetPCC(i,j);//need to change
				os<<m_vProtein_name[i]<<'\t'<<m_vProtein_name[j]<<'\t'<<GetECC(i,j)<<'\t'<<GetPCC(i,j)<<'\t'<<m_vvWeight[i][j]<<endl;
			}
		}
		os.close();
}

double CPPI::Similarity( const Complex &c1,const Complex &c2 )
{
	int p1=0,p2=0,same=0;
	int sz1=c1.size(),sz2=c2.size();
	assert(sz1>0&&sz2>0);
	//since c1 & c2 are sorted,we use O(n) method to calculate intersecton set
	while(p1<sz1&&p2<sz2)
	{
		if(c1[p1]==c2[p2])
		{
			++same;
			++p1;
			++p2;
		}
		else if(c1[p1]<c2[p2])
			++p1;
		else
			++p2;
	}

	return (double)same*same/(sz1*sz2);
}

bool CPPI::Test( const Complex &c, vector<Complex> &v )
{
	if(c.size()<=1)
		return 0;
	for(int i=0;i<v.size();++i)
		if(Similarity(v[i],c)>m_argR)
			return 0;
	return 1;
}

double CPPI::Similarity2( const Complex &c1,const Complex &c2 )
{
	int p1=0,p2=0,same=0;
	int sz1=c1.size(),sz2=c2.size();
	assert(sz1>0&&sz2>0);
	//since c1 & c2 are sorted,we use O(n) method to calculate intersecton set

	while(p1<sz1&&p2<sz2)
	{
		if(c1[p1]==c2[p2])
		{
			++same;
			++p1;
			++p2;
		}
		else if(c1[p1]<c2[p2])
			++p1;
		else
			++p2;
	}
	//两个集合完全相同或其中一个完全包含在另一个集合内
	return same==min(c1.size(),c2.size());
}

bool CPPI::Test2( const Complex &c, vector<Complex> &v )
{
	if(c.size()<=1)
		return 0;
	for(int i=0;i<v.size();++i)
		if(Similarity2(v[i],c))
			return 0;
	return 1;
}

void CPPI::FindCoreEss()
{
	//获取核
	for(int i=0;i<m_num;++i)
	{
		if(m_is_ess[i])
		{
			Complex cp;
			cp.push_back(i);
			for(int id:m_vvAdjcency_list[i])
			{
				if(m_is_ess[id]||m_vvWeight[id][i]>m_argT)
					cp.push_back(id);
			}
			sort(cp.begin(),cp.end());
			m_M1.push_back(cp);
		}
	}
}

void CPPI::ExpandEss()
{
	for(int i=0;i<m_M1.size();++i)
	{
		Complex tmp(m_M1[i]);
		vector <bool> vis(m_num);
		for(int id:tmp)
			vis[id]=1;
		for(int n1:tmp)
		{
			for(int n2:m_vvAdjcency_list[n1])
			{
				if(!vis[n2]&&m_vvWeight[n1][n2]>m_argT)
				{
					vis[n2]=1;
					m_M1[i].push_back(n2);
				}
			}
		}
		sort(m_M1[i].begin(),m_M1[i].end());
		assert(unique(m_M1[i].begin(),m_M1[i].end())==m_M1[i].end());
	}

	sort(m_M1.begin(),m_M1.end(),
		[](const Complex &c1,const Complex &c2){return c1.size()>c2.size();});
	vector <Complex> vtmp;
	for(int i=0;i<m_M1.size();++i)
	{
		//满足两个条件才加入复合物集合
		//test：与前面所有复合物的相似性小于R
		//test2：不完全包含在另一个集合内
		if(Test(m_M1[i],vtmp)&&Test2(m_M1[i],vtmp))
			vtmp.push_back(m_M1[i]);
	}
	m_M1=vtmp;
}

int CPPI::TwoHopDegree( int id )
{
	int ans=m_vvAdjcency_list[id].size();
	for(int n1:m_vvAdjcency_list[id])
		ans+=m_vvAdjcency_list[n1].size();
	return ans;
}

double CPPI::Density( const Complex &c )
{
	int sz=c.size();
	if(sz<=1) return 0;
	int e=0;
	for(int i=0;i<sz;++i)
		for(int j=i+1;j<sz;++j)
			if(m_vvAdjcency_matrix[c[i]][c[j]])
				++e;
	return 2*e/(sz*(sz-1));
}

void CPPI::FindCoreNoEss()
{
	vector <bool> used(m_num);
	for(auto& comp:m_M1)
		for(int id:comp)
			used[id]=1;

	vector <int> H;
	for(int i=0;i<m_num;++i)
	{
		if(!used[i])
			H.push_back(i);
	}

	m_vtwohop.resize(m_num);
	for(int i=0;i<m_num;++i)
		m_vtwohop[i]=TwoHopDegree(i);

	sort(H.begin(),H.end(),
		[this](const int &a,const int &b){return m_vtwohop[a] > m_vtwohop[b];});  //H按二级度排序
	for(int id:H)
		if(!used[id])
		{
			Complex tmp;
			tmp.push_back(id);
			vector <int> neighbors;
			for(int n1:m_vvAdjcency_list[id])
				if(!used[n1])
					neighbors.push_back(n1);
			sort(neighbors.begin(),neighbors.end(),
				[this](const int &a,const int &b){return m_vvAdjcency_list[a].size()>m_vvAdjcency_list[b].size();}); //邻居按度排序
			for(int n1:neighbors)
			{
				tmp.push_back(n1);
				if(Density(tmp)<0.7)    //加入新点后使密度<0.7则剔除
					tmp.pop_back();
			}
			for(int id:tmp)
				used[id]=1;
			if(tmp.size()>2)
				m_M2.push_back(tmp);
		}
}



void CPPI::ExpandNoEss()
{
	for(int i=0;i<m_M2.size();++i)
	{
		Complex tmp(m_M2[i]);
		vector <bool> vis(m_num);
		for(int id:tmp)
			vis[id]=1;
		for(int n1:tmp)
		{
			for(int n2:m_vvAdjcency_list[n1])
			{
				if(!vis[n2])  //遍历复合物的邻居，且如果其未被访问过
				{
					vis[n2]=1;
					int deg=0;
					for(int k:tmp)   //遍历复合物内的点，计算内度
					{
						if(m_vvAdjcency_matrix[n2][k])
							++deg;
					}
					//如果有一半以上的点相连，则加入复合物核
					if(2*deg>tmp.size())
						m_M2[i].push_back(n2);
				}
			}
		}
		sort(m_M2[i].begin(),m_M2[i].end());
		assert(unique(m_M2[i].begin(),m_M2[i].end())==m_M2[i].end());
	}
}
