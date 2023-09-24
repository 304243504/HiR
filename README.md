# HiR analysis pipeline

## Introduction
This is a specialized pipeline under development to analysis RNA-RNA interaction information captured by an explicit bridge linker. The pipeline consists of three main steps:
1. linker identification
2. sequence alignment
3. Cluster Identification

## Usage
Python3 environment is a must.
Before excuting the HiR pipeline, you need to create hisat2 genome index and onfigure environment variables of bedtools.
In addation, create the genome_info file using the following command line:

`python getGeneElement.py -i your_specie_annotation.gtf -o genome.element.bed`


To excute the HiR pipeline, simply run:

`sh Run.sh ${prefix} linker_file Int[length(linker)*0.9]`


## Main references
1. Cai, Z., Cao, C., Ji, L., Ye, R., Wang, D., Xia, C., Wang, S., Du, Z., Hu, N., Yu, X., et al. (2020). RIC-seq for global in situ profiling of RNA-RNA spatial interactions. Nature 582, 432–437. https://doi.org/10.1038/s41586-0202249-1.
2. Xiao, Q., Huang, X., Zhang, Y., Xu, W., Yang, Y., Zhang, Q., Hu, Z., Xing, F., Sun, Q., Li, G., & Li, X. (2022). The landscape of promoter-centred RNA-DNA interactions in rice. Nature plants, 8(2), 157–170. https://doi.org/10.1038/s41477-021-01089-4.




