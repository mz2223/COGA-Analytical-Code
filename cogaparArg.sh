#!/bin/bash 

#SBATCH -A r01370
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --time=3:59:00
#SBATCH --mem=60gb
#SBATCH --mail-user=mz22@iu.edu
#SBATCH -J cogapar
#SBATCH -o cogaparArg.o
#SBATCH -e cogaparArg.e

date
module load python
python cogaparArg.py 1 22
tr -d ' ' < allout1-22.csv > trimmed.csv
mv trimmed.csv "ea-$(basename "$PWD").csv"
#tar -czf "ea-$(basename "$PWD").csv.tar.gz" trimmed.csv
date
