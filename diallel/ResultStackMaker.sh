#!/bin/bash
source PostRunConfig.sh
START=$(`echo date`)

##---------------------------------------------
# Make stack file and lists
##---------------------------------------------
mkdir "${DIR}/output/"
touch "${DIR}/output/stackedTreatDiMIMPs.csv"
touch "${DIR}/output/stackedPostTreatDiMIMPs.csv"
touch "${DIR}/output/stackedPostPredTreatDiMIMPs.csv"
touch "${DIR}/MIMP_listCompleted.out"
touch "${DIR}/MIMP_listPostCompleted.out"
touch "${DIR}/MIMP_listPostPredCompleted.out"

##---------------------------------------------
# Get list of existing MIMP files
##---------------------------------------------
find "$DIR" -type f -name 'StackedTreatDi.csv' >> "${DIR}/MIMP_listCompleted.out"
NUMFILES="$(cat ${DIR}/MIMP_listCompleted.out | sed '/^\s*$/d' | wc -l)"
echo "The number of MIMP files is: ${NUMFILES}"

##---------------------------------------------
# Get list of existing MIMP Posterior files
##---------------------------------------------
find "$DIR" -type f -name 'StackedTreatDiPost.csv' >> "${DIR}/MIMP_listPostCompleted.out"
NUMFILES="$(cat ${DIR}/MIMP_listPostCompleted.out | sed '/^\s*$/d' | wc -l)"
echo "The number of MIMP posterior files is: ${NUMFILES}"

##-----------------------------------------------------
# Get list of existing MIMP Posterior Predictive files
##-----------------------------------------------------
find "$DIR" -type f -name 'StackedTreatDiPostPred.csv' >> "${DIR}/MIMP_listPostPredCompleted.out"
NUMFILES="$(cat ${DIR}/MIMP_listPostPredCompleted.out | sed '/^\s*$/d' | wc -l)"
echo "The number of MIMP posterior predictive files is: ${NUMFILES}"

##---------------------------------------------
# Make header; loop through files
##---------------------------------------------
FIRSTFILE=$(head -n 1 "${DIR}/MIMP_listCompleted.out")
head -n 1 $FIRSTFILE > "${DIR}/output/stackedTreatDiMIMPs.csv"

while read line; do
    ./seqlines.pl 2 END $THIN $line >> "${DIR}/output/stackedTreatDiMIMPs.csv"
    echo "stacked ${line}"
done < "${DIR}/MIMP_listCompleted.out"

echo "Result thinning and stacking completed."

##---------------------------------------------
## Make header; loop through posterior files
##---------------------------------------------
FIRSTFILE=$(head -n 1 "${DIR}/MIMP_listPostCompleted.out")
head -n 1 $FIRSTFILE > "${DIR}/output/stackedPostTreatDiMIMPs.csv"

while read line; do
    ./seqlines.pl 2 END $THIN $line >> "${DIR}/output/stackedPostTreatDiMIMPs.csv"
    echo "stacked ${line}"
done < "${DIR}/MIMP_listPostCompleted.out"

echo "Posterior thinning and stacking completed."
END=$(`echo date`)

##---------------------------------------------
## Make header; loop through posterior files
##---------------------------------------------
FIRSTFILE=$(head -n 1 "${DIR}/MIMP_listPostPredCompleted.out")
head -n 1 $FIRSTFILE > "${DIR}/output/stackedPostPredTreatDiMIMPs.csv"

while read line; do
    ./seqlines.pl 2 END $THIN $line >> "${DIR}/output/stackedPostPredTreatDiMIMPs.csv"
    echo "stacked ${line}"
done < "${DIR}/MIMP_listPostPredCompleted.out"

echo "Posterior predictive thinning and stacking completed."
END=$(`echo date`)

##---------------------------------------------
echo "STACKING START : $START"
echo "STACKING END   : $END"
