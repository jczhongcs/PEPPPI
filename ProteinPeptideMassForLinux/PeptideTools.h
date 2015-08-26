#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include <map>
#include <string>
using namespace std;
/*

template <typename T, typename U>
class create_map
{
private:
	std::map<T, U> m_map;
public:
	create_map(const T& key, const U& val)
	{
		m_map[key] = val;
	}

	create_map<T, U>& operator()(const T& key, const U& val)
	{
		m_map[key] = val;
		return *this;
	}

	operator std::map<T, U>()
	{
		return m_map;
	}
};*/
class CPeptideTools {
public:
    CPeptideTools();

	static map<char,double> mapMonoisotopic;

	static double dHydroneMass;



	static double massQuantity(const string source)
	{
		double dmassQuantify=0.0;
		for (int i = 0;i<source.length();i++)
		{
			dmassQuantify+=mapMonoisotopic[source.at(i)];
		}

		return dmassQuantify;

	}
	//Ëã·¨,calculate edit distance
	static int ldistance(const string source,const string target)
	{
		//step 1

		int n=source.length();
		int m=target.length();
		if (m==0) return n;
		if (n==0) return m;
		//Construct a matrix
		typedef vector< vector<int> >  Tmatrix;
		Tmatrix matrix(n+1);
		for(int i=0; i<=n; i++)  matrix[i].resize(m+1);

		//step 2 Initialize

		for(int i=1;i<=n;i++) matrix[i][0]=i;
		for(int i=1;i<=m;i++) matrix[0][i]=i;

		//step 3
		for(int i=1;i<=n;i++)
		{
			const char si=source[i-1];
			//step 4
			for(int j=1;j<=m;j++)
			{

				const char dj=target[j-1];
				//step 5
				int cost;
				if(si==dj){
					cost=0;
				}
				else{
					cost=1;
				}
				//step 6
				const int above=matrix[i-1][j]+1;
				const int left=matrix[i][j-1]+1;
				const int diag=matrix[i-1][j-1]+cost;
				matrix[i][j]=min(above,min(left,diag));

			}
		}//step7
		return matrix[n][m];
	}

	static void getnextval(const char* T,int next[])
	{
		int j=0,k=-1;
		next[0]=-1;
		while(T[j]!='\0')
		{
			if(k==-1||T[j]==T[k])
			{
				j++; k++;
				if(T[j]!=T[k]) next[j]=k;
				else next[j]=next[k];
			}
			else k=next[k];
		}
	}

	//ÉèÔÚ×Ö·û´®SÖÐ²éÕÒÄ£Ê½´®T,ÈôS[m]!=T[n],È¡T[n]Ä£Ê½º¯ÊýÖµnext[n],
	//Èç¹ûnext[n]=-1,±íÊ¾S[m]ºÍT[0]¼ä½Ó±È½Ï¹ýÁË,²»ÏàµÈ,ÏÂÒ»´Î±È½Ï S[m+1] ºÍT[0]
	//Èç¹ûnext[n]=0 ±íÊ¾±È½Ï¹ý³ÌÖÐ²úÉúÁË²»ÏàµÈ£¬ÏÂÒ»´Î±È½Ï S[m] ºÍT[0]
	//Èç¹ûnext[n]= k >0 µ«k<n, ±íÊ¾,S[m]µÄÇ°k¸ö×Ö·ûÓëTÖÐµÄ¿ªÊ¼k¸ö×Ö·ûÒÑ¾­¼ä½Ó±È½ÏÏàµÈÁË£¬ÏÂÒ»´Î±È½ÏS[m]ºÍT[k]ÏàµÈÂð£¿
	//////////////////////////////////////////////////////////////////////////
	//parameterÔ­Ê¼´®£¬Ä£Ê½´®
	static int KMP(char *Text,const char *Pattern)
	{
		if(!Text||!Pattern||Text=='\0'||Pattern=='\0') return -1;
		int len=strlen(Pattern);
		int *next=new int[len+1];
		getnextval(Pattern,next);
		int i=0,j=0,index=0;
		while(Text[i]!='\0'&&Pattern[j]!='\0')
		{
			if(Text[i]==Pattern[j])
			{
				i++; j++;
			}
			else
			{
				index+=j-next[j];
				if(next[j]!=-1) j=next[j];
				else
				{
					i++; j=0;
				}
			}
		}
		if(Pattern[j]=='\0') return index;
		return -1;
	}

	//////////////////////////////////////////////////////////////////////////
	//the efficience is better than kmp
	static int SUNDAYFindSigle(char *text, char *patt){ 
		size_t  temp[256]; 
		size_t  *shift = temp; 
		size_t  i, patt_size = strlen(patt), text_size = strlen(text); 
//		cout << "size : " << patt_size << endl; 
		for( i=0; i < 256; i++ ) *(shift+i) = patt_size+1; 
		for( i=0; i < patt_size; i++ )  
			*(shift + (unsigned char)(*(patt+i))) = patt_size-i; 
		//shift['s']=6?,shitf['e']=5 ???? 
		size_t  limit = text_size-patt_size+1; 
		for( i=0; i < limit; i += shift[ text[i+patt_size] ] ) 
			if( text[i] == *patt ){ 
				char *match_text = text+i+1; 
				size_t  match_size = 1; 
				do{// Êä³öËùÓÐÆ¥ÅäµÄÎ»ÖÃ 
					if( match_size == patt_size ) return i;//cout << "the 	NO. is " << i << endl; 
				}while( (*match_text++) == 	patt[match_size++] ); 
			}

	    return -1;//no find pattern string
//			cout << endl; 
	}


	//////////////////////////////////////////////////////////////////////////
	//the efficience is better than kmp
	static int SUNDAYFind(char *text, char *patt, int iLine,int iIsform,ofstream& out1,string& name){ 
		size_t  temp[256]; 
		size_t  *shift = temp; 
		size_t  i, patt_size = strlen(patt), text_size = strlen(text); 
		//		cout << "size : " << patt_size << endl; 
		for( i=0; i < 256; i++ ) *(shift+i) = patt_size+1; 
		for( i=0; i < patt_size; i++ )  
			*(shift + (unsigned char)(*(patt+i))) = patt_size-i; 
		//shift['s']=6?,shitf['e']=5 ???? 
		size_t  limit = text_size-patt_size+1; 
		for( i=0; i < limit; i += shift[ text[i+patt_size] ] ) 
			if( text[i] == *patt ){ 
				char *match_text = text+i+1; 
				size_t  match_size = 1; 
				do{// Êä³öËùÓÐÆ¥ÅäµÄÎ»ÖÃ 
//					if( match_size == patt_size ) return i;//cout << "the 	NO. is " << i << endl; 
					if( match_size == patt_size ) out1 << iLine<<'\t'<< name<<'\t'<<iIsform<<'\t'<<text<<'\t'<<patt<<'\t'<<i<<endl;; 

				}while( (*match_text++) == 	patt[match_size++] ); 
			}

			return -1;//no find pattern string
			//			cout << endl; 
	}


	/*int main(){
		string s;
		string d;
		cout<<"source=";
		cin>>s;
		cout<<"diag=";
		cin>>d;
		int dist=ldistance(s,d);
		cout<<"dist="<<dist<<endl;
	}}*/
  
};