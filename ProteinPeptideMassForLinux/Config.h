#pragma once
#include <string>
#include <fstream>
using namespace std;
class Config
{
public:
	Config(void);
	~Config(void);


	static string strPPIFileName;
	static string strSequenceFileName;
	static string strGoldensetFileName;
	static string strKeyProteinFileName;

	static string strPeptideProphetResult;
	static char cSplitChar;

	static int iRunTimes;
	
	static string strOutResultFileName;


	static void ReadPara()
	{
		ifstream config_file("para.ini");
	
		string strFromPara;
		while(config_file>>strFromPara)
		{
			if(strFromPara=="strPPIFileName")
				config_file>>strPPIFileName;
			else if(strFromPara=="strSequenceFileName")
				config_file>>strSequenceFileName;
			else if(strFromPara=="strGoldensetFileName")
				config_file>>strGoldensetFileName;
			else if(strFromPara=="strKeyProteinFileName")
				config_file>>strKeyProteinFileName;
			else if(strFromPara=="strPeptideProphetResult")
				config_file>>strPeptideProphetResult;
			else if(strFromPara=="cSplitChar")
				config_file>>cSplitChar;
			else if(strFromPara=="iRunTimes")
				config_file>>iRunTimes;
			else if(strFromPara=="strOutResultFileName")
				config_file>>strOutResultFileName;
			


		}

	}
};

