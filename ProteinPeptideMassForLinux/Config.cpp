#include "stdafx.h"
#include "Config.h"




string Config::strPPIFileName("C:\\Users\\Administrator\\Desktop\\massInference\\20140711result\\data04\\DIP20150429Filter.txt");
string Config::strSequenceFileName("j:\\mass\\mass_yeast\\sc_SGD_0604.fasta");
string Config::strGoldensetFileName("C:\\Users\\Administrator\\Desktop\\massInference\\20140711result\\data02GoldenSet.txt");
string Config::strPeptideProphetResult("E:\\workP1\\ProteinLP\\real_data\\yeastData04XTgt05_6.txt");
char   Config::cSplitChar(',');
int Config::iRunTimes=500;
string Config::strOutResultFileName("outFirstResultYeast1");

Config::Config(void)
{
}


Config::~Config(void)
{
}
