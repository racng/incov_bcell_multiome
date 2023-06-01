#!/bin/bash

DATADIR=/ssd1/rng/data/incov_bcell_multiome/
REFDIR=/home/rng/reference/
BCDIR=/home/rng/proj/incov_bcell_multiome/output/cell_calling/
VCFDIR=/ssd1/rng/output/demultiplex/merge-wgs
POOLS='R5 R6 R7 R8'

THD='32'

for ID in $POOLS; 
do
    # Navigate to library folder
    cd $DATADIR/$ID/
    mkdir vireo
    ## Vireo using genotypes.  Two samples pooled.    
    vireo -t GT -N 2  \
    --vartrixData=vartrix/vt_alt.mtx,vartrix/vt_ref.mtx,${BCDIR}/accepted_barcodes_outer_${ID}.csv,${VCFDIR}/${ID}.var.exonp250_hg38.vcf.gz  \
    -d ${VCFDIR}/${ID}.var.exonp250_hg38.vcf.gz -o vireo 
done