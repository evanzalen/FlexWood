#!/bin/bash -l

set -ex

# environment variables
proj=u2021023

mail="elena.vanzalen@umu.se"

in=$(realpath ../)
out=$(realpath ../MultiQC)

if [ ! -d $out ]; then
	mkdir -p $out
fi

module load bioinfo-tools MultiQC

sbatch -A $proj --mail-user=$mail -o $in/MultiQC/multiqc.out -e $in/MultiQC/multiqc.err \
-J $(basename $in) ../UPSCb-common/pipeline/runMultiQC.sh $in $out