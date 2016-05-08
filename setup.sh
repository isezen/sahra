#!/usr/bin/env bash

R --slave<<EOF
if (!"rwind" %in% rownames(installed.packages()))
  devtools::install_github("isezen/rwind")
if (!"rsezen" %in% rownames(installed.packages()))
  devtools::install_github("isezen/rsezen")
source("code/nchelper.r")
download_ncep_R2()
nc2rds_all()
EOF
