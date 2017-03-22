#!/bin/bash

date

Rscript PostRunTemp.R --args --configfile=settings.config.D1 &
Rscript PostRunTemp.R --args --configfile=settings.config.D2 &
Rscript PostRunTemp.R --args --configfile=settings.config.D3 &
Rscript PostRunTemp.R --args --configfile=settings.config.D4 &

echo "============================================================"
echo "Submitted PostRun!"
echo "============================================================"

date


