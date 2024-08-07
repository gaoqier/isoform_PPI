# Set locale to C
Sys.setlocale("LC_ALL", "C")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

packages <- c("tidyverse", "broom", "dbplyr", "ggplot2", "modelr", "ragg", "reprex", "rvest", "stringr", "tidyr")

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Check and install missing packages
lapply(packages, install_if_missing)

# Load the packages
lapply(packages, library, character.only = TRUE)

# List of packages needed
packages <- c("tidyverse", "parallel")

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Check and install missing packages
lapply(packages, install_if_missing)

# Load the packages
lapply(packages, library, character.only = TRUE)



library(tidyverse)
library(parallel)

samples="SRR3649792"
sj.list = mclapply(samples, function(sample){   
     sj_sp = data.table::fread(paste0("./work/star/",sample,"/",sample,".SJ.out.tab"), data.table=F)
     if(nrow(sj_sp)==0) return(NULL)
     sj_sp = sj_sp %>% 
       dplyr::mutate(`coord.intron`=paste0(V1,":",V2,":",V3)) %>% 
       dplyr::select(`coord.intron`, V7)
     colnames(sj_sp)[2] = sample
     return(sj_sp)
},mc.cores=10)
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
