#!/bin/bash
code_dir="new_sparsity_test_barycentric"

for submission_dir in $code_dir/*/*/*/*/*/; do
#    cd $submission_dir
    scancel $(find $submission_dir/. -type f -iname "slurm*" | grep -o -P '(?<=slurm-).*(?=.out)')
#    pwd
#    cd -
done

