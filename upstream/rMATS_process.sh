id=SRR3649792
ref_idx_star=../basic/gtf/star_index
ref_gtf=../basic/gtf/gencode.v40.annotation.gtf
cd /gaoqier/MARVEL/proj1
mkdir -p /gaoqier/MARVEL/proj1/work/rmats/${id}
echo ${id} is running......
echo "./work/star/${id}/${id}.Aligned.sortedByCoord.out.bam" > ./work/star/${id}/BAM_fls.txt

rmats.py \
--b1 ./work/star/${id}//BAM_fls.txt \
--gtf ${ref_gtf} \
--od ./work/rmats/${id} \
--tmp ./work/rmats/${id}/tmp \
-t paired \
--readLength 100 \
--variable-read-length \
--nthread 12 \
--statoff
# readLength参考如下
# cat ./work/fq/${id}_1.fastq | awk '{if(NR%4==2) print NR"\t"$0"\t"length($0)}'
