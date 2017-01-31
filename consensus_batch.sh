#!/bin/bash
# author: Astride Audirac
# date: 01/06/15
# script to run a batch 
# 

INPUTDIR=$1;
SEED=$2;
OUTPUTDIR=$3;
TMPDIR=$4;
PROC=$5;
$j=0;

echo sh consensus_batch.sh  $SEED $OUTPUTDIR $TMPDIR $PROC

#mkdir $OUTPUTDIR
cd $INPUTDIR
for i in $(ls);
	do
		$j=$j+1;
		cat "i"|grep -B1 $SEED> "'i'.${j}.fasta"; 
	done
	


