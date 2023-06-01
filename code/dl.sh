#!/bin/bash
aws s3 cp  s3://isb-covid-single-cell/processed_data/cellranger_arc_out/ ./ --recursive --exclude "*" --include "R*/outs/per_barcode_metrics.csv"
aws s3 cp  s3://isb-covid-single-cell/processed_data/cellranger_arc_out/ ./ --recursive --exclude "*" --include "R[5-8]/outs/atac_possorted_bam.bam"
aws s3 cp  s3://isb-covid-single-cell/processed_data/cellranger_arc_out/ ./ --recursive --exclude "*" --include "R[5-8]/outs/atac_possorted_bam.bam.bai"
aws s3 cp  s3://isb-covid-single-cell/processed_data/cellranger_arc_out/ ./ --recursive --exclude "*" --include "R[5-8]/outs/gex_possorted_bam.bam"
aws s3 cp  s3://isb-covid-single-cell/processed_data/cellranger_arc_out/ ./ --recursive --exclude "*" --include "R[5-8]/outs/gex_possorted_bam.bam.bai"
