# incov_bcell_multiome

## Code


|File|Description|
|-|-|
|dl.sh|Script to download barcode metrics and ATAC and GEX BAM files.|
|run_script.sh bcell|Script to run [snakemake workflow](https://github.com/racng/snakemake-merge-wgs) to generate pooled VCF fileGenerate a VCF file per sample|
|config_bcell.yaml|Config file for snakemake workflow|
|merge_bam.sh|Merge RNA and ATAC BAM files.|
|config.py|Plotting settings and paths.|
|call_cells.py|Perform joint calling of cells based on UMI count and ATAC peak counts.|
|vartrix.sh|Extract single cell variants from merged BAM file using [VarTrix](https://github.com/10XGenomics/vartrix) |
|vireo.sh|Run [vireo](https://github.com/single-cell-genetics/vireo) demultiplexing using vartrix output and sample VCF file.|
|check_demult.py|Compare donor ratios with experimental ratios.|

