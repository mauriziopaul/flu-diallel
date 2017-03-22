#!/bin/bash

date

Rscript PostPredImpute.R --args --configfile=settings.config.D3-Diplo6 &&
Rscript ImputeMatch.R --args --configfile=settings.config.D3-Diplo6 &&
Rscript PileMake.R --args --configfile=settings.config.D3-Diplo6 &&
bash BatchSend.sh

echo "============================================================"
echo "Submitted PostPredImpute, PileMake and BatchSend!"
echo "Please call MIMQ_PostRun.sh to run ResultStacker and PostRun ONLY when piles are complete."
echo "============================================================"

date
