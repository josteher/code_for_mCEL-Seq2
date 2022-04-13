#! /bin/bash

g=$(find $(pwd) -name \*_run_mapping.sh)

for i in $g;do echo $i; sbatch $i; sleep 1; done
