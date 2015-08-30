# PEPPPI
Protein inference from the integration of tandem MS data and PPI networks
#*PEPPPI*

**Written by** Jiancheng Zhong (jczhong@csu.edu.cn)  
[Jianxin Wang Lab, Central South University](http://netlab.csu.edu.cn/)

**Please cite:**

---

**Current version:** 0.1.0

Support for Windows

##Summary
In PPI networks, interacting (thus co-existing) proteins for a specific function are often neighbors. We present our proposed method, called *PEPPPI*, by introducing the extended bipartite model combined with PPI networks. The PEPPPI iteratively computes the probability of a shared peptide belonging to a specific protein. Compared with the original bipartite model, our PEPPPI integrates PPI networks with tandem MS data to infer proteins. 


##Installation

#### Requirements
visual studio 2012 C++

##### Install
```
git clone git@github.com:jczhongcs/PEPPPI.git
cd PEPPPI
open ProteinPeptideMassForLinux.sln with visual studio 2012
run *PEPPPI* that was built with Visual Studio 2012
```
##Usage
User can modify the para.ini to set the parameters for the PEPPPI, which include the PPIFileName, SequenceFileName, GoldenSetFileName(optional),PeptideProphetResult, Splitchar, RunTimes, and OutResultFileName.
For example(yeast datasets)
strPPIFileName	../RealDatasets/yeast/DIP20131031.txt
strSequenceFileName	../RealDatasets/yeast/sc_SGD_0604.fasta
strGoldensetFileName	../RealDatasets/yeast/GoldenSet.txt
strPeptideProphetResult	../RealDatasets/yeast/yeastData04XTgt05_6.txt
cSplitChar	,
iRunTimes	500
strOutResultFileName	../RealDatasets/yeast/outFirstResultYeast1
