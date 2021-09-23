#$ -S /bin/bash
#$ -N Filt1Dead
#$ -V
#$ -pe smp 10
#$ -cwd
#$ -l mem_free=64G
#$ -o 06_Filt1Dead_STDOUT.log
#$ -e 06_Filt1Dead_STDERR.log

/vcu_gpfs2/home/mccbnfolab/harrell-pdx-scrnaseq/Harrell_SingleCellSequencing/bash/03_filterCells.sh -f /vcu_gpfs2/home/mccbnfolab/harrell-pdx-scrnaseq/config/06_filterDeadCells_singleSamples_grch38_070121.list -d /vcu_gpfs2/home/mccbnfolab/harrell-pdx-scrnaseq/raw/human_only_deadfilt/
