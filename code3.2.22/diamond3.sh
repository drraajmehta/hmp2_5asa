#!/bin/bash
#SBATCH -n 4                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 6-23:30          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=32000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid

source /n/huttenhower_lab/tools/hutlab/src/hutlabrc.sh
hutlab load centos7/python3/humann3/3.0.0-alpha-devel

diamond blastp -q candidate5asa.fa -d /net/fs-huttenhower01/srv/export/huttenhower_lab/share/data/humann2_databases/version_3.0.0-alpha/uniref90/uniref90_201901.dmnd -o out4.tsv --max-target-seqs 0 --evalue 10 --id 25 --sensitive


#wget https://www.uniprot.org/uniprot/Q00267.fasta
#wget https://www.uniprot.org/uniprot/P77567.fasta
#cat /n/home02/rmehta/*.fasta > candidate5asa.fa




#good explanation of query coverage https://www.biostars.org/p/342445/

#diamond/0.9.5-fasrc01