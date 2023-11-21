#$ -S /bin/bash
#$ -N cutadapt
#$ -V
#$ -pe smp 20
#$ -cwd
#$ -l mem_free=10G
#$ -o 01_cutadapt_STDOUT.log
#$ -e 01_cutadapt_STDERR.log

cd /vcu_gpfs2/home/harrell_lab/bulkRNASeq/

INFILE=''

mkdir -p 01_cutadapt

cat $INFILE | while read s1 s2
do
	base1=`basename $s1 .fastq.gz`
	base2=`basename $s2 .fastq.gz`
	cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -q 15 -o 01_cutadapt/$base1.trimmed.fastq.gz -p 01_cutadapt/$base2.trimmed.fastq.gz $s1 $s2
done
