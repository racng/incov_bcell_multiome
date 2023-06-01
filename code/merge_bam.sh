#!/bin/bash
# conda activate vireo

DATADIR=/ssd1/rng/data/incov_bcell_multiome/
POOLS='R5 R6 R7 R8'

THD='32'

for ID in $POOLS; do
	cd $DATADIR/$ID/outs/
	
	## Merge RNA and ATAC bam files.
	samtools merge merged.bam -@ $THD gex_possorted_bam.bam atac_possorted_bam.bam
	samtools index merged.bam 
done
