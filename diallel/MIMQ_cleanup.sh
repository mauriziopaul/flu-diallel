#!/bin/bash
find . -type f -name "TreatDi_MIMP.RData" -delete &
find . -type f -name "PartSaveListTBSR5.RData" -delete &
find . -type f -name "SaveTBSR5.RData" -delete &
find . -type f -name "*PostPredImpute.RData" -delete &
find . -type f -name "StackedTreatDi*.csv*" -delete &
find . -type f -name "TreatDi*.csv" -delete &
find . -type f -name "AutoDiallelChain*.bin" -delete &
find . -type f -name "Pile*.sh" -delete &
find . -type f -name "BatchSend.sh" -delete &
find . -type f -name "MIM*_list*.out" -delete &
find . -type f -name "PostRunConfig.sh" -delete &
find . -type f -name "output_postrunner.out" -delete &
find . -type f -name "*PostPredImpute.RData" -delete &
find . -type f -name "ImputeCensorStack.csv" -delete &
find . -type f -name "InnerRun_Completed.txt" -delete &

#find . -type f -name "MIM*_List.txt" -delete &
#find . -type f -name "stackedTreatDiMIMPs.csv" -delete &
#find . -type f -name "stackedPostTreatDiMIMPs.csv" -delete &
#find . -type f -name "stackedPostPredTreatDiMIMPs.csv" -delete &
#find . -type f -name "new.mip.table.csv" -delete &
#find . -type f -name "PSqTableStacked.csv" -delete &
