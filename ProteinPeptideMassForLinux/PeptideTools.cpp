#include "stdafx.h"
#include "PeptideTools.h"
//myClass.cpp
map<char,double>  CPeptideTools::mapMonoisotopic = map<char,double>();
double CPeptideTools::dHydroneMass=19.01783;


CPeptideTools::CPeptideTools()
{

	mapMonoisotopic.insert(make_pair('A',71.03711));
	mapMonoisotopic.insert(make_pair('R',156.10111));
	mapMonoisotopic.insert(make_pair('N',114.04293));
	mapMonoisotopic.insert(make_pair('D',115.02694));
	mapMonoisotopic.insert(make_pair('C',103.00919));
	mapMonoisotopic.insert(make_pair('E',129.04259));
	mapMonoisotopic.insert(make_pair('Q',128.05858));
	mapMonoisotopic.insert(make_pair('G',57.02146));
	mapMonoisotopic.insert(make_pair('H',137.05891));
	mapMonoisotopic.insert(make_pair('I',113.08406));
	mapMonoisotopic.insert(make_pair('L',113.08406));
	mapMonoisotopic.insert(make_pair('K',128.09496));
	mapMonoisotopic.insert(make_pair('M',131.04049));
	mapMonoisotopic.insert(make_pair('F',147.06841));
	mapMonoisotopic.insert(make_pair('P',97.05276));
	mapMonoisotopic.insert(make_pair('S',87.03203));
	mapMonoisotopic.insert(make_pair('T',101.04768));
	mapMonoisotopic.insert(make_pair('W',186.07931));
	mapMonoisotopic.insert(make_pair('Y',163.06333));
	mapMonoisotopic.insert(make_pair('V',99.06841));
	mapMonoisotopic.insert(make_pair('m',15.9949));
	mapMonoisotopic.insert(make_pair('c',57.02165));
}