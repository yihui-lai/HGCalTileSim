#*******************************************************************************
 #
 #  Filename    : test_qeury.py
 #  Description : Unit testing for query testings
 #  Author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
 #
#*******************************************************************************
import UserUtils.EDMUtils.EDMFileList as edmlist

print edmlist.ListEDMDataSets('SingleMuon','*Run2016*','MINIAOD')
print edmlist.ListFileInDataset('/SingleMuon/Run2016B-PromptReco-v1/MINIAOD')
