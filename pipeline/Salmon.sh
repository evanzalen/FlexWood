#!/bin/bash -l

set -ex

# environment variables
proj=u2021023

mail="elena.vanzalen@umu.se"

in=$(realpath ../data/raw)
out=$(realpath ../)
ref=/mnt/picea/storage/reference/Populus-tremula/v2.2/indices/salmon1.2.1


start=7
end=7
#mem=256 
export PATH=${PATH}:$(realpath ../UPSCb-common/kogia/scripts)

module load bioinfo-tools SortMeRNA trimmomatic
for f in `find $in -name "*_1.fq.gz"`; do
bash $(realpath ../UPSCb-common/pipeline/runRNASeqPreprocessing.sh) -s $start -e $end -S $ref $proj $mail $f ${f/_1.fq.gz/_2.fq.gz} $out
done
