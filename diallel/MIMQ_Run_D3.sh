#!/bin/bash

date

Rscript PostPredImpute.R --args --configfile=settings.config.D3 &&
Rscript ImputeMatch.R --args --configfile=settings.config.D3 &&
Rscript PileMake.R --args --configfile=settings.config.D3 &&
bash BatchSend.sh

echo "============================================================"
echo "Submitted PostPredImpute, PileMake and BatchSend!"
echo "Please call MIMQ_PostRun.sh to run ResultStacker and PostRun ONLY when piles are complete."
echo "============================================================"

date
