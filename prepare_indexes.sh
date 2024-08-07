#!/bin/bash

# 定义文件路径
GENOME_FA="GRCh38.primary_assembly.genome.fa"
GTF_FILE="gencode.v40.annotation.gtf"
STAR_INDEX_DIR="star_index"
RSEM_INDEX_DIR="rsem_index/hg38"

# 生成 STAR 基因组索引
STAR --runMode genomeGenerate \
     --genomeDir ${STAR_INDEX_DIR} \
     --genomeFastaFiles ${GENOME_FA} \
     --sjdbGTFfile ${GTF_FILE} \
     --runThreadN 12

# 生成 RSEM 转录本索引
rsem-prepare-reference \
     --gtf ${GTF_FILE} \
     --star -p 12 \
     ${GENOME_FA} \
     ${RSEM_INDEX_DIR}
