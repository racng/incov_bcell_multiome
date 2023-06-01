#!/bin/bash
# conda activate vireo

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
    # Make folder for vartrix outputs
    mkdir vartrix

    /home/rng/prog/vartrix_linux --threads $THD --scoring-method coverage  \
        -v ${VCFDIR}/${ID}.var.exonp250_hg38.vcf.gz  \
        -b outs/merged.bam  \
        -f ${REFDIR}/hg38.fa  \
        -c ${BCDIR}/accepted_barcodes_outer_${ID}.csv  \
        -o vartrix/vt_alt.mtx --ref-matrix vartrix/vt_ref.mtx \
        --out-variants vartrix/vt_var.mtx 
done
