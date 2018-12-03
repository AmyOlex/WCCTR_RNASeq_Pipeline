#!/bin/bash
#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M alolex@vcu.edu
#PBS -N splitByOrg
#PBS -j oe
#PBS -o splitByOrg.log


INPUT=/home/alolex/src/ChuckHarrell_BrainMetastasis_12-2016/RNAseq_Processing_Pipeline/file_lists/input05_toSplitByOrg.list
OUTDIR=/home/sequencing/data/WorkData/Harrell_AllBatches_OfficialWorkflow/05_splitByOrg

mkdir -p $OUTDIR

echo "Starting Splitting "`date`

## Human lines are identified by "ENST" and Mouse lines are identified by "ENSMUST".

HUMAN="ENST"
MOUSE="ENSMUST"

cat $INPUT | while read file
do
	base=`basename $file .quant.sf`
	head -n 1 $file > head.tmp
	grep $HUMAN $file > human.tmp
	grep $MOUSE $file > mouse.tmp
	cat head.tmp human.tmp > $OUTDIR/$base.human.quant.sf
	cat head.tmp mouse.tmp > $OUTDIR/$base.mouse.quant.sf

	rm head.tmp
	rm human.tmp
	rm mouse.tmp
done

echo "Finished Splitting "`date`
