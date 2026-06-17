#!/bin/bash
code_dir="voltage_1.00eV"

for submission_dir in $code_dir/*/*/*/*/; do
#    cd $submission_dir
     rm $(find $submission_dir/. -type f -iname "*.dat")
     rm $(find $submission_dir/. -type f -iname "slurm*")
     rm $(find $submission_dir/. -type f -iname "*.p")
     rm $(find $submission_dir/. -type f -iname "*.npy")
#    pwd
#    cd -
done
