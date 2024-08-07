#2.1 创建marvel对象
library(tidyverse)
library(MARVEL)
library(readxl)
# (1) 细胞表型分组

df.pheno <- readxl::read_excel("./SJ_phenoData.xlsx")

# (2) Gene文件
df.tpm <- read.table("./work/merge/rsem_tpm.txt", sep="\t", header=TRUE)
dim(df.tpm)


# (3) GTF文件
gtf <- data.table::fread("../basic/gtf/gencode.v40.annotation.gtf", sep="\t", header=FALSE, na.strings="NA",quote="\"") %>% 
  as.data.frame()

# (4) Gene类型
gene.feature = subset(gtf, V3=="gene") 
df.tpm.feature <- str_split(gene.feature$V9, ";", simplify=T) %>% 
  as.data.frame() %>%
  dplyr::mutate(gene_id=str_match(V1, '"(.*)"')[,2]) %>% 
  dplyr::mutate(gene_short_name=str_match(V3, '"(.*)"')[,2]) %>%      
  dplyr::mutate(gene_type=str_match(V2, '"(.*)"')[,2]) %>% 
  dplyr::select(gene_id, gene_short_name, gene_type) %>% 
  dplyr::filter(gene_id %in% df.tpm$gene_id)
df.tpm.feature = df.tpm.feature[match(df.tpm$gene_id, df.tpm.feature$gene_id),]
rownames(df.tpm.feature) = seq(nrow(df.tpm.feature))

# (5) SJ文件
sj = data.table::fread("./work/merge/SJ.txt", sep="\t", header=TRUE, na.strings="NA") %>% 
  as.data.frame()
dim(sj)

# (6) AS文件
df.feature.list = lapply(c("SE", "MXE", "RI", "A5SS", "A3SS"), function(x){
  # x = "SE"
  df.feature = read.table(paste0("./work/merge/as_",x,"_featureData.txt"), 
                          sep="\t", header=TRUE, na.strings="NA") %>% 
    dplyr::distinct(tran_id, .keep_all=TRUE) %>% 
    dplyr::left_join(df.tpm.feature)
  return(df.feature)
})
names(df.feature.list) <- c("SE", "MXE", "RI", "A5SS", "A3SS")

# (7) Intron文件
df.intron.counts <- data.table::fread("./work/merge/intron_count_by_region.txt", 
                                      sep="\t", header=TRUE, na.strings="NA") %>% 
  as.data.frame()
dim(df.intron.counts)
## MARVEL object
marvel <- CreateMarvelObject(SpliceJunction=sj,
                             SplicePheno=df.pheno,
                             SpliceFeature=df.feature.list,
                             IntronCounts=df.intron.counts,
                             GeneFeature=df.tpm.feature,
                             Exp=df.tpm,
                             GTF=gtf
)

#2.2 预处理步骤
marvel <- TransformExpValues(MarvelObject=marvel,
                             offset=1,
                             transformation="log2",
                             threshold.lower=1
)
# 如果基因表达水平整体较低，则不考虑该基因的可变剪切事件
# min.expr:判定单个细胞是否表达基因的阈值
# min.cells: 基于上述条件，判定细胞群表达该基因的阈值
marvel <- DetectEvents(MarvelObject=marvel,
                       min.cells=5, #细胞群的表达百分比
                       min.expr=0.5,   #认为细胞表达该基因的阈值
                       track.progress=FALSE,
                       EventType="AFE"
)
marvel <- DetectEvents(MarvelObject=marvel,
                       min.cells=5,
                       min.expr=0.5,
                       track.progress=FALSE,
                       EventType="ALE"
)
length(marvel$SpliceFeature)
# 7
# CoverageThreshold: 支持该剪切事件的最小reads数，小于该阈值标记为NA
# UnevenCoverageMultiplier: 针对SE与MXE
marvel <- CheckAlignment(MarvelObject=marvel, level="SJ")
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     UnevenCoverageMultiplier=10,
                     EventType="SE"
)

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     UnevenCoverageMultiplier=10,
                     EventType="MXE"
)

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="RI",
                     thread=4  # only support RI
)

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="A5SS"
)

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="A3SS"
)

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="AFE"
)

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="ALE"
)
marvel <- CheckAlignment(MarvelObject=marvel, level="SJ")
marvel <- CheckAlignment(MarvelObject=marvel, level="splicing")
marvel <- CheckAlignment(MarvelObject=marvel, level="gene")
marvel <- CheckAlignment(MarvelObject=marvel, level="splicing and gene")

#根据情况加
#index.1 <- which(df.pheno$cell.type %in% c("iPSC", "Endoderm"))
#index.2 <- which(df.pheno$qc.seq=="pass")

#index <- intersect(index.1, index.2)
#sample.ids <- df.pheno[index, "sample.id"]
#marvel <- SubsetSamples(MarvelObject=marvel,
#                        sample.ids=sample.ids
#)

save(marvel, file="./MARVEL2.RData")

# 后续以官方提供的示例数据学习
# load("../Data/MARVEL.RData")





#2.3 细胞群AS特征分析
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno

# Define sample ids
sample.ids <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]

# Tabulate expressed events
marvel <- CountEvents(MarvelObject=marvel,
                      sample.ids=sample.ids,
                      min.cells=25
)

# Output (1): Plot
marvel$N.Events$Plot
