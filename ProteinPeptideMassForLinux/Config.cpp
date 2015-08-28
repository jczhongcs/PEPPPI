#include "stdafx.h"
#include "Config.h"




string Config::strPPIFileName("../RealDatasets/yeast/DIP20131031.txt");
string Config::strSequenceFileName("../RealDatasets/yeast/sc_SGD_0604.fasta");
string Config::strGoldensetFileName("../RealDatasets/yeast/GoldenSet.txt");
string Config::strPeptideProphetResult("../RealDatasets/yeast/yeastData04XTgt05_6.txt");
char   Config::cSplitChar(',');
int Config::iRunTimes=500;
string Config::strOutResultFileName("../RealDatasets/yeast/outFirstResultYeast1");

Config::Config(void)
{
}


Config::~Config(void)
{
}
