cd /gaoqier/MARVEL/proj1
id=SRR3649792
mkdir -p ./work/rsem/${id}
rsem-calculate-expression \
--bam \
-p 12 \
./work/star/${id}/${id}.Aligned.toTranscriptome.out.bam \
../basic/gtf/rsem_index/hg38 \
./work/rsem/${id}/${id}