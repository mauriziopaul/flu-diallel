#!/bin/bash

date
Rscript PostRunTemp.R --args --configfile=settings.config.D1 &
Rscript PostRunTemp.R --args --configfile=settings.config.D2 &
Rscript PostRunTemp.R --args --configfile=settings.config.D3 &
Rscript PostRunTemp.R --args --configfile=settings.config.D4 &
Rscript PostRunTemp.R --args --configfile=settings.config.D1-Diplo6 &
Rscript PostRunTemp.R --args --configfile=settings.config.D2-Diplo6 &
Rscript PostRunTemp.R --args --configfile=settings.config.D3-Diplo6 &
Rscript PostRunTemp.R --args --configfile=settings.config.D4-Diplo6 &
date


