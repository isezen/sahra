#!/usr/bin/env bash

R --slave<<EOF
 options(cl.cores = 7) # number of cores
 source("code/calcor.r")
EOF
