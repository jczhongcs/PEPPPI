#include "stdafx.h"
#include "CalProteinProbLogisticRegression.h"





CCalProteinProbLogisticRegression::~CCalProteinProbLogisticRegression(void)
{
}

void CCalProteinProbLogisticRegression::calculateProteinProbability()
{
	double dProbability=1.0;
	for (CMassPeptideXcorr* pPeptideXcorr:m_pProtein->m_vSupportPeptide)//Loop SupportPeptide
	{
		int iPeptideSharedProtein=1;
		double dOriginPeptideSP=1.0;
		if (pPeptideXcorr->m_bUnique==false)//shared peptide (sequest data)
		{
			double dOtherProteinPeptideValue=0.0;

			
			dOriginPeptideSP=getSameMassSPValue( pPeptideXcorr);//get the peptide sigma sp value, mass match to peptide 

			vector<double> vdPeptideSP(pPeptideXcorr->m_vProtein.size(),0.0);//save the sibling peptide value
			for (int j=0;j<pPeptideXcorr->m_vProtein.size();j++)//Loop share peptide to Protein
			{
				vector<CMassPeptideXcorr*> vMassPeptide = m_pListProtein->getSupportPeptideFromProteinStr(pPeptideXcorr->m_vProtein[j]);
				double dPeptideSP=0.0;
				for (int j1=0;j1<vMassPeptide.size();j1++)//every protein has some peptides 
				{
					if (vMassPeptide[j1]->m_strPeptideSequence!=pPeptideXcorr->m_strPeptideSequence)//except for the original peptideSequence
					{
						dPeptideSP += getSameMassSPValue( vMassPeptide[j1]);
					}

				}
				dPeptideSP/=(pPeptideXcorr->m_vProtein.size()-1);

				vdPeptideSP.push_back(dPeptideSP);
				if (pPeptideXcorr->m_vProtein[j]!=m_pProtein->m_strIPIName)
				{
					dOtherProteinPeptideValue+=dPeptideSP;
				}

			}
			dOriginPeptideSP=abs(dOriginPeptideSP-dOtherProteinPeptideValue)/dOriginPeptideSP;
		

		}

		dProbability=dProbability*(1-dOriginPeptideSP);//no exist probability

	}
	m_pProtein->m_dProbability=pow((1-dProbability),-log(m_pProtein->m_vSupportPeptide.size()/m_pProtein->m_ivalidTheoryPeptideCount));
}
