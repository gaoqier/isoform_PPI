#安装包
install.packages("ggplot2")
install.packages("Matrix")
install.packages("plyr")
install.packages("scales")
install.packages("devtools")
library(devtools)
install_github("wenweixiong/MARVEL")
library(MARVEL)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")
BiocManager::install("Biostrings")
BiocManager::install("BSgenome")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("clusterProfiler")
BiocManager::install("GenomicRanges")
BiocManager::install("IRanges")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("Rtools")
install.packages("assertthat")


# 设置镜像源
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# 增加下载超时时间
options(timeout = max(1200, getOption("timeout")))
# 再次尝试安装包
BiocManager::install("phastCons100way.UCSC.hg38")

install.packages("wiggleplotr_1.18.0.tar.gz", repos=NULL, type="source")

# Load MARVEL package
library(MARVEL)

# Load adjunct packages for selected MARVEL features
# General data processing, plotting
library(ggnewscale)
library(ggrepel)
library(parallel)
library(reshape2)
library(stringr)
install.packages("textclean")
library(textclean)

# Dimension reduction analysis
install.packages("factoextra")
library(factoextra)
library(FactoMineR)

# Modality analysis
install.packages("fitdistrplus")
library(fitdistrplus)

# Differential splicing analysis
library(kSamples)
library(twosamples)

# Gene ontology analysis
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# Nonsense-mediated decay (NMD) analysis
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.NCBI.GRCh38)

# Load adjunct packages for this tutorial
library(data.table)
library(ggplot2)
library(gridExtra)

#Sample metadata
path <- "Data/SJ/"
file <- "SJ_phenoData.txt"
df.pheno <- read.table(paste(path, file, sep=""), header=TRUE, stringsAsFactors=FALSE, na.strings="NA")

head(df.pheno)

