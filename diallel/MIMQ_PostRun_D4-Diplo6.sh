#!/bin/bash

date

bash ResultStackMakerGZ.sh &&
source PostRunConfig.sh && Rscript PostRun.R --args --configfile=settings.config.D4-Diplo6

echo "============================================================"
echo "Submitted StackCaller and PostRun!"
echo "============================================================"

date


