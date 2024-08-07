library(tidyverse)
library(parallel)
install.packages("readxl")
library(readxl)
# 读取 XLSX 文件
sample_meta <- readxl::read_excel("./SJ_phenoData.xlsx")
samples <- sample_meta$sample.id




#step1: 汇总剪切位点矩阵
library(tidyverse)
library(parallel)
library(data.table)
library(dplyr)
samples="SRR3649792"
cl <- makeCluster(detectCores() - 1)  # 使用所有可用内核减去一个
sj.list <- parLapply(cl, samples, function(sample) {
  # 在并行计算节点上加载所需包
  library(data.table)
  library(dplyr)
  
  sj_sp <- data.table::fread(paste0("./work/star/", sample, "/", sample, ".SJ.out.tab"), data.table = FALSE)
  if (nrow(sj_sp) == 0) return(NULL)
  sj_sp <- sj_sp %>%
    dplyr::mutate(`coord.intron` = paste0(V1, ":", V2, ":", V3)) %>%
    dplyr::select(`coord.intron`, V7)
  colnames(sj_sp)[2] <- sample
  return(sj_sp)
})
names(sj.list) = samples
sj.list = sj.list[!unlist(lapply(sj.list, is.null))]

sj_res = sj.list[[1]]
for(x in sj.list[-1]){
  sj_res = dplyr::full_join(sj_res, x)
}

for(sample in setdiff(samples,names(sj.list))){
  sj_res[,sample]=NA
}
dim(sj_res)

write.table(sj_res, file=paste0("./work/merge/SJ.txt"),
            quote=F, row.names=F, sep="\t")

#step2:整理AS注释结果
# https://github.com/lishensuo/utils/blob/main/marvel/rscript_Anno_SJ_rMATs.R
source(paste0("../basic/rscript_Anno_SJ_rMATs.R"))

# 定义 AS 类型
AS_types = c("SE", "RI", "MXE", "A3SS", "A5SS")

# 并行处理每种 AS 类型
for(AS_type in AS_types) {
  print(AS_type)
  
  fls_AS = list.files(paste0("./work/rmats"), pattern=paste0("fromGTF.", AS_type, ".txt"),
                      recursive = TRUE, full.names = TRUE)
  
  if (AS_type == "SE") {
    func = se_func
  } else if (AS_type == "RI") {
    func = ri_func 
  } else if (AS_type == "MXE") {
    func = mxe_func 
  } else if (AS_type == "A3SS") {
    func = a3ss_func 
  } else if (AS_type == "A5SS") {
    func = a5ss_func
  }
  
  cl <- makeCluster(detectCores() - 1)
  
  # 确保每个节点加载所需的包
  clusterExport(cl, c("func", "fls_AS", "se_func", "ri_func", "mxe_func", "a3ss_func", "a5ss_func"))
  clusterEvalQ(cl, {
    library(data.table)
    library(dplyr)
  })
  
  df_AS <- parLapply(cl, fls_AS, function(fls) {
    rmats_AS <- data.table::fread(fls, data.table = FALSE)
    rmats_AS <- rmats_AS %>%
      dplyr::mutate(tran_id = func(.)) %>%
      dplyr::rename(gene_id = GeneID, gene_short_name = geneSymbol) %>%
      dplyr::select(tran_id, gene_id, gene_short_name)
    return(rmats_AS)
  }) %>% do.call(rbind, .) %>% 
    dplyr::distinct()
  
  stopCluster(cl)
  
  write.table(df_AS, row.names = FALSE, sep = "\t", quote = FALSE,
              file = paste0("./work/merge/as_", AS_type, "_featureData.txt"))
}





#step3：汇总intron矩阵
library(data.table)
library(dplyr)
library(parallel)
# 创建并行集群
cl <- makeCluster(detectCores() - 1)

# 确保每个节点加载所需的包和导出样本列表
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
})
clusterExport(cl, "samples")

int_mt.list = parLapply(cl, samples, function(sample){
  # sample = samples[1]
  int_mt = data.table::fread(paste0("./work/intron/",sample,"/intron_count.txt")) 
  int_mt = int_mt %>% 
    dplyr::mutate(`coord.intron`=paste(V1,V2,V3, sep=":")) %>% 
    dplyr::group_by(coord.intron) %>% 
    dplyr::summarize(count=sum(V5)) %>% 
    dplyr::select(coord.intron, count) %>% as.data.frame()
  colnames(int_mt)[2] = sample
  return(int_mt)
})

intron_res = int_mt.list[[1]]
for(x in int_mt.list[-1]){
  intron_res = dplyr::full_join(intron_res, x)
}
dim(intron_res)

write.table(intron_res, file=paste0("./work/merge/intron_count_by_region.txt"),
            quote=F, row.names=F, sep="\t")


#汇总gene矩阵
library(data.table)
library(dplyr)
library(parallel)
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
})
clusterExport(cl, "samples")

gene_mt.list = parLapply(cl, samples, function(sample){
  fls=paste0("./work/rsem/",sample,"/",sample,".genes.results")
  gene_mt = data.table::fread(fls, data.table=F) %>% 
    dplyr::select(gene_id, TPM)
  colnames(gene_mt)[2] = sample
  return(gene_mt)
})
gene_res = gene_mt.list[[1]]
for(x in gene_mt.list[-1]){
  gene_res = dplyr::full_join(gene_res, x)
}
dim(gene_res)
write.table(gene_res, file=paste0("./work/merge/rsem_tpm.txt"),,
            quote=F, row.names=F, sep="\t")
