#include "stdafx.h"
#include "EnzymeProteinPeptide.h"
#include <string>
#include <boost/regex.hpp>

CEnzymeProteinPeptide::CEnzymeProteinPeptide(void)
{
	//////////////////////////////////////////////////////////////////////////
	  /****************************************
		 *                                       *
		 *            Chop'N'Spice               *
		 *                                       *
		 *       Main Configuration File         *
		 *                                       *
		 *****************************************/
	mapEnzymeRegex["Arg-C proteinase"]  = "(.*?R|(?<=R).+?$)";
	mapEnzymeRegex["Asp-N endopeptidase"]  = "(.+?(?=D)|D.*?$)";
	mapEnzymeRegex["Asp-N endopeptidase (+N-terminal Glu)"]  = "(.+?(?=D|E)|[DE].*?$)";
	mapEnzymeRegex["Chymotrypsin - high specificity"]  = "(.*?([FY](?!P)|W(?!M|P))|(?<=F|Y).+?$)";
	mapEnzymeRegex["Chymotrypsin - low specificity"]  = "(.*?([FLY](?!P)|W(?!M|P)|M(?!P|Y)|H(?!D|M|P|W))|(?<=F|L|Y).+?$)";
	mapEnzymeRegex["CNBr - high excess"]  = "(.*?M|(?<=M).+?$)";
	mapEnzymeRegex["CNBr - low excess"]  = "(.*?M(?!S|T)|(?<=M)[^ST].*?$)";
	mapEnzymeRegex["Enterokinase"]  = "(.*?[DN]{3}K|(?<=[DN]{3}).+?$)";
	mapEnzymeRegex["Factor Xa"]  = "(.*?[AFGILTVM]{1}[DE]{1}GR|(?<=[AFGILTVM]{1}[DE]{1}GR).+?$)";
	mapEnzymeRegex["Formic acid"]  = "(.*?D|(?<=D).+?$)";
	mapEnzymeRegex["Glu-C endopeptidase"]  = "(.*?E(?!P)|(?<=E)[^P].*?$)";
	mapEnzymeRegex["Glu-C endopeptidase (N-term C or D)"]  = "(.*?[DE](?!P)|(?<=D|E)[^P].*?$)";
	mapEnzymeRegex["Granzyme B"]  = "(.*?IEPD|(?<=IEPD).+?$)";
	mapEnzymeRegex["Hydroxylamine"]  = "(.*?N(?=G)|(?<=N)G.*?$)";
	mapEnzymeRegex["Iodosobenzoic acid"]  = "(.*?W|(?<=W).+?$)";
	mapEnzymeRegex["Lys-C proteinase"]  = "(.*?K(?!P)|(?<=K)[^P].*?$)";
	mapEnzymeRegex["Lys-C proteinase (+C-term P)"]  = "(.*?K|(?<=K).+?$)";
	mapEnzymeRegex["NTCB - 2-nitro-5-thiocyanobenzoic acid"]  = "(.+?(?=C)|C.*?$)";
	mapEnzymeRegex["Pepsin - pH 1.3"]  = "(.*?[^HKR][^P]([^R](?=[FLWY][^P])|[FLWY](?=.[^P]))|(?<=[^HKR][^P][^R])[FLWY][^P].*?$|(?<=[^HKR][^P][FLWY]).[^P].+?$)";
	mapEnzymeRegex["Pepsin - pH > 2"]  = "(.*?[^HKR][^P]([^R](?=[FL][^P])|[FL](?=.[^P]))|(?<=[^HKR][^P][^R])[FL][^P].*?$|(?<=[^HKR][^P][FL]).[^P].+?$)";
	mapEnzymeRegex["Proteinase K"]  = "(.*?[AEFILTVWY]|(?<=[AEFILTVWY]).+?$)";
	mapEnzymeRegex["TEV protease"]  = "(.*?ENLYFQ(?=G)|(?<=ENLYFQ)G.*?$)";
	mapEnzymeRegex["Thrombin"]  = "(.*?(GR(?=G)|[AFGILTVM][AFGILTVWA]PR(?=[^DE]{2}))|(?<=GR)G.*?$|(?<=[AFGILTVM][AFGILTVWA]PR)[^DE]{2}.*?$)";
	mapEnzymeRegex["Thermolysin"]  = "(.*?[^DE](?=A|F|I|L|M|V)|(?<=[^DE])[AFILMV].*?$)";
	mapEnzymeRegex["Trypsin - N-term K or R"]  = "(.*?[K,R]|(?<=[K,R]).+?$)";
	mapEnzymeRegex["Trypsin - N-term K or R, not C-term P"]  = "(.*?[K,R](?!P)|(?<=[K,R])[^P].*?$)";
	mapEnzymeRegex["Trypsin - advanced model"]  = "(.*?([^CD]K(?!P)|WK(?=P)|DK(?!D)|CK(?!D|H|Y)|[^CR]R(?!P)|MR(?=P)|CR(?!K)|RR(?!H|R))|(?<=K|R).+?$)";
	mapEnzymeRegex["none"]  = "(^.*$)";

}


CEnzymeProteinPeptide::~CEnzymeProteinPeptide(void)
{
}



void CEnzymeProteinPeptide::SplitProteinWithEnzyme( string& strProteinName,string& strProteinSequence,string& strEnzymeType,map<string,string>& m_mapPeptideProtein )
{
	//set EnzymeRegex

	const boost::regex pattern(mapEnzymeRegex[strEnzymeType], boost::regex::perl|boost::regex::icase);

	boost::sregex_iterator it(strProteinSequence.begin(),strProteinSequence.end(),pattern);

	boost::sregex_iterator end;

	while (it!=end) 
	{
		m_mapPeptideProtein[(*it)[1]] = strProteinName;

		++it;

	}
}

void CEnzymeProteinPeptide::SplitProteinWithEnzyme( string& strProteinSequence,string& strEnzymeType,vector<string>& vPeptides,int& ivalidPeptideCount )
{
	//set EnzymeRegex

	const boost::regex pattern(mapEnzymeRegex[strEnzymeType], boost::regex::perl|boost::regex::icase);

	boost::sregex_iterator it(strProteinSequence.begin(),strProteinSequence.end(),pattern);

	boost::sregex_iterator end;

	int iBz=0;
	while (it!=end) 
	{

		iBz++;
		vPeptides.push_back((*it)[1]);

		//////////////////////////////////////////////////////////////////////////
		//这里可以考虑剪切长度数小于N的情况
		if ((*it)[1].length()>5)//小于5的peptide不是有效peptide
		{
			ivalidPeptideCount++;
		}

		++it;

	}

	if (iBz==0)//若没有剪切
	{
		vPeptides.push_back(strProteinSequence);
	}
	if (ivalidPeptideCount==0)
	{
		ivalidPeptideCount=1;
	}
}
