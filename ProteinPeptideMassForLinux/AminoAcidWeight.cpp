#include "stdafx.h"
#include "AminoAcidWeight.h"

/*

	'A' => '71.037114',
	'C' => '103.009185',
	'D' => '115.026943',
	'E' => '129.042593',
	'F' => '147.068414',
	'G' => '57.021464',
	'H' => '137.058912',
	'I' => '113.084064',
	'J' => '0.0',
	'K' => '128.094963',
	'L' => '113.084064',
	'M' => '131.040485',
	'N' => '114.042927',
	'P' => '97.052764',
	'Q' => '128.058578',
	'R' => '156.101111',
	'S' => '87.032028',
	'T' => '101.047679',
	'V' => '99.068414',
	'W' => '186.079313',
	'Y' => '163.063329',
	'HYDROGEN' => '1.007825035',
	'CARBON' => '12.0',
	'NITROGEN' => '14.003074',
	'OXYGEN' => '15.99491463',
	'ELECTRON' => '0.000549',
	'PROTON' => '1.007276035',
	'H2O' => '18.0105647'
*/


const double CAminoAcidWeight::m_rc[26]={ /* A */  0.8, /* B */  0.0, /* C */ -0.8, /* D */ -0.5, /* E */  0.0, /* F */ 10.5,
	/* G */ -0.9, /* H */ -1.3, /* I */  8.4, /* J */  0.0, /* K */ -1.9, /* L */  9.6,
	/* M */  5.8, /* N */ -1.2, /* O */  0.0, /* P */  0.2, /* Q */ -0.9, /* R */ -1.3,
	/* S */ -0.8, /* T */  0.4, /* U */  0.0, /* V */  5.0, /* W */ 11.0, /* X */  0.0,
	/* Y */  4.0, /* Z */  0.0};


const double CAminoAcidWeight::m_rcnt[26]={  /* A */ -1.5, /* B */  0.0, /* C */  4.0, /* D */  9.0, /* E */  7.0, /* F */ -7.0,
	/* G */  5.0, /* H */  4.0, /* I */ -8.0, /* J */  0.0, /* K */  4.6, /* L */ -9.0,
	/* M */ -5.5, /* N */  5.0, /* O */  0.0, /* P */  4.0, /* Q */  1.0, /* R */  8.0,
	/* S */  5.0, /* T */  5.0, /* U */  0.0, /* V */ -5.5, /* W */ -4.0, /* X */  0.0,
	/* Y */ -3.0, /* Z */  0.0};

const double CAminoAcidWeight::m_nt[3] = {0.42, 0.22, .05};

const double CAminoAcidWeight::m_minMass = 1.0;

const double CAminoAcidWeight::m_maxMass = 640000.0;


CAminoAcidWeight::CAminoAcidWeight()
{

	//////////////////////////////////////////////////////////////////////////
	//Initailize MonoisotopicWeight
/*
	m_mapMonoisotopicWeight["A"] = 71.037114;
	m_mapMonoisotopicWeight["C"] = 103.009185;
	m_mapMonoisotopicWeight["D"] = 115.026943;
	m_mapMonoisotopicWeight["E"] = 129.042593;
	m_mapMonoisotopicWeight["F"] = 147.068414;
	m_mapMonoisotopicWeight["G"] = 57.021464;
	m_mapMonoisotopicWeight["H"] = 137.058912;
	m_mapMonoisotopicWeight["I"] = 113.084064;
	m_mapMonoisotopicWeight["J"] = 0.0;
	m_mapMonoisotopicWeight["K"] = 128.094963;
	m_mapMonoisotopicWeight["L"] = 113.084064;
	m_mapMonoisotopicWeight["M"] = 131.040485;
	m_mapMonoisotopicWeight["N"] = 114.042927;
	m_mapMonoisotopicWeight["P"] = 97.052764;
	m_mapMonoisotopicWeight["Q"] = 128.058578;
	m_mapMonoisotopicWeight["R"] = 156.101111;
	m_mapMonoisotopicWeight["S"] = 87.032028;
	m_mapMonoisotopicWeight["T"] = 101.047679;
	m_mapMonoisotopicWeight["V"] = 99.068414;
	m_mapMonoisotopicWeight["W"] = 186.079313;
	m_mapMonoisotopicWeight["Y"] = 163.063329;*/



	m_mapMonoisotopicWeight["G"]= 57.0214636;
	m_mapMonoisotopicWeight["A"]= 71.0371136;
	m_mapMonoisotopicWeight["S"]= 87.0320282;
	m_mapMonoisotopicWeight["P"]= 97.0527636;
	m_mapMonoisotopicWeight["V"]= 99.0684136;
	m_mapMonoisotopicWeight["T"]=101.0476782;
	m_mapMonoisotopicWeight["C"]=103.0091854 +57.021464; //Carboxamidomethyl (on Cysteine); Added by Zengyou (28-02-2008)
	m_mapMonoisotopicWeight["L"]=113.0840636;
	m_mapMonoisotopicWeight["I"]=113.0840636;
	m_mapMonoisotopicWeight["X"]=113.0840636;
	m_mapMonoisotopicWeight["N"]=114.0429272;
	m_mapMonoisotopicWeight["O"]=114.0793126;
	m_mapMonoisotopicWeight["B"]=114.5349350;
	m_mapMonoisotopicWeight["D"]=115.0269428;
	m_mapMonoisotopicWeight["Q"]=128.0585772;
	m_mapMonoisotopicWeight["K"]=128.0949626;
	m_mapMonoisotopicWeight["Z"]=128.5505850;
	m_mapMonoisotopicWeight["E"]=129.0425928;
	m_mapMonoisotopicWeight["M"]=131.0404854;
	m_mapMonoisotopicWeight["H"]=137.0589116;
	m_mapMonoisotopicWeight["F"]=147.0684136;
	m_mapMonoisotopicWeight["R"]=156.1011106;
	m_mapMonoisotopicWeight["Y"]=163.0633282;
	m_mapMonoisotopicWeight["W"]=186.0793126;

	m_mapMonoisotopicWeight["HYDROGEN"] = 1.007825035;
	m_mapMonoisotopicWeight["CARBON"] = 12.0;
	m_mapMonoisotopicWeight["NITROGEN"] = 14.003074;
	m_mapMonoisotopicWeight["OXYGEN"] = 15.99491463;
	m_mapMonoisotopicWeight["ELECTRON"] = 0.0005485;
	m_mapMonoisotopicWeight["PROTON"] = 1.007276035;
	m_mapMonoisotopicWeight["H2O"] = 18.0105647;
	m_mapMonoisotopicWeight["PHOSPORUS"]= 30.9737633;
	m_mapMonoisotopicWeight["SULPHUR"]= 31.9720718;


	//Initailize AverageWeight

/*
	m_mapAverageWeight["A"] = 71.0779;
	m_mapAverageWeight["C"] = 103.1429;
	m_mapAverageWeight["D"] = 115.0874;
	m_mapAverageWeight["E"] = 129.114;
	m_mapAverageWeight["F"] = 147.1739;
	m_mapAverageWeight["G"] = 57.0513;
	m_mapAverageWeight["H"] = 137.1393;
	m_mapAverageWeight["I"] = 113.1576;
	m_mapAverageWeight["J"] = 0.0;
	m_mapAverageWeight["K"] = 128.1723;
	m_mapAverageWeight["L"] = 113.1576;
	m_mapAverageWeight["M"] = 131.1961;
	m_mapAverageWeight["N"] = 114.1026;
	m_mapAverageWeight["P"] = 97.1152;
	m_mapAverageWeight["Q"] = 128.1292;
	m_mapAverageWeight["R"] = 156.1857;
	m_mapAverageWeight["S"] = 87.0773;
	m_mapAverageWeight["T"] = 101.1039;
	m_mapAverageWeight["V"] = 99.1311;
	m_mapAverageWeight["W"] = 186.2099;
	m_mapAverageWeight["Y"] = 163.1733;


	m_mapAverageWeight["h"]=  1.00794;  / * hydrogen * /
	m_mapAverageWeight["o"]= 15.9994;   / * oxygen * /
	m_mapAverageWeight["c"]= 12.0107;   / * carbon * /
	m_mapAverageWeight["n"]= 14.00674;  / * nitrogen * /
	m_mapAverageWeight["p"]= 30.973761; / * phosporus * /
	m_mapAverageWeight["s"]= 32.066;    / * sulphur * /*/

	m_mapAverageWeight["G"]= 57.05192;
	m_mapAverageWeight["A"]= 71.07880;
	m_mapAverageWeight["S"]= 87.07820;
	m_mapAverageWeight["P"]= 97.11668;
	m_mapAverageWeight["V"]= 99.13256;
	m_mapAverageWeight["T"]=101.10508;
	m_mapAverageWeight["C"]=103.13880 +57.0513; /* 103.1448, 103.14080 */ //Carboxamidomethyl(on Cysteine);Added by Zengyou(28-02-2008)
	m_mapAverageWeight["L"]=113.15944;
	m_mapAverageWeight["I"]=113.15944;
	m_mapAverageWeight["X"]=113.15944;
	m_mapAverageWeight["N"]=114.10384;
	m_mapAverageWeight["O"]=114.14720;
	m_mapAverageWeight["B"]=114.59622;
	m_mapAverageWeight["D"]=115.08860;
	m_mapAverageWeight["Q"]=128.13072;
	m_mapAverageWeight["K"]=128.17408;
	m_mapAverageWeight["Z"]=128.62310;
	m_mapAverageWeight["E"]=129.11548;
	m_mapAverageWeight["M"]=131.19256; /* 131.19456 131.1986 */
	m_mapAverageWeight["H"]=137.14108;
	m_mapAverageWeight["F"]=147.17656;
	m_mapAverageWeight["R"]=156.18748;
	m_mapAverageWeight["Y"]=163.17596;
	m_mapAverageWeight["W"]=186.21320;

	m_mapAverageWeight["HYDROGEN"] = 1.00794;
	m_mapAverageWeight["CARBON"] = 12.0107;
	m_mapAverageWeight["NITROGEN"] = 14.00674;
	m_mapAverageWeight["OXYGEN"] = 15.9994;
	m_mapAverageWeight["ELECTRON"] = 0.0005485;
	m_mapAverageWeight["PROTON"] = 1.007391;
	m_mapAverageWeight["H2O"] = 18.01528;
	m_mapAverageWeight["PHOSPORUS"]= 30.973761; /* phosporus */
	m_mapAverageWeight["SULPHUR"]= 32.066;    /* sulphur */




}


CAminoAcidWeight::~CAminoAcidWeight()
{
}

double CAminoAcidWeight::getMonoisotopicWeight( string& strPeptideSequence )
{
	double dWeight=0.0;
	for (int i=0;i<strPeptideSequence.size();i++)
	{
		dWeight+=m_mapMonoisotopicWeight[strPeptideSequence.substr(i,1)];
	}
	//need add other modification
	return dWeight;
}

double CAminoAcidWeight::getAverageWeight( string& strPeptideSequence )
{
	double dWeight=0.0;
	for (int i=0;i<strPeptideSequence.size();i++)
	{
		dWeight+=m_mapAverageWeight[strPeptideSequence.substr(i,1)];
	}
	//need add other modification
	return dWeight;
}
//////////////////////////////////////////////////////////////////////////
//ÊèË®ÐÔ
//source: ScoreRegularization
double CAminoAcidWeight::hydrophobicity( string& bytes, int start, int n )
{
	double kl=1;
	if (n<10)
	{
		kl=1-0.027*(10-n);
	}else
	{
		if (n>20)
		{
			kl=1/(1+0.015*(n-20));//revision from paper's author
			//kl=1-0.014*(n-20)//as published in paper
		}
	}
	double h=0;
	for (int i=0;i<n;i++)
	{
		char c=(char)bytes[i+start];
		h += m_rc[c - 'A'];
		if(i<3)
			h+=m_nt[i]*m_rcnt[c-'A'];
	}
	h*=kl;
	if(h<38)
		return h ;
	else
		return h-0.3*(h-38);

}
