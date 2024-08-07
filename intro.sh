cd /gaoqier/MARVEL/proj1
id=SRR3649792
ref_idx_star=../basic/gtf/star_index
ref_gtf=../basic/gtf/gencode.v40.annotation.gtf
mkdir -p ./work/intron/${id}
echo ${id} is running......

samtools view -H ./work/star/${id}/${id}.Aligned.sortedByCoord.out.bam | \
  grep SQ | cut -f 2 | awk '{ sub(/^SN:/, ""); print;}' > ./work/star/${id}/sorted_chr_in_bam.txt
# 预定义R脚本，生成bedtools计算所需的两个文件：染色体大小以及内含子坐标
# https://github.com/lishensuo/utils/blob/main/marvel/rscript_bedtools_input.R
Rscript ../basic/rscript_bedtools_input.R ${id}


bedtools coverage \
-g ./work/intron/${id}/hg38.chrom.sizes.txt \
-split \
-sorted \
-a ./work/intron/${id}/RI_Coordinates_sorted.bed \
-b ./work/star/${id}/${id}.Aligned.sortedByCoord.out.bam > \
./work/intron/${id}/intron_count.txt \
-d
