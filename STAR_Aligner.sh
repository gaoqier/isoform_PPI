cd /gaoqier/MARVEL/proj1
id=SRR3649792
mkdir -p ./work/star/${id}
ref_idx_star=../basic/gtf/star_index
ref_gtf=../basic/gtf/gencode.v40.annotation.gtf
echo ${id} is runing......
#第一轮
STAR --runThreadN 12 \
--genomeDir ${ref_idx_star} \
--readFilesIn ./work/fq/${id}/${id}.fastq \
--outFileNamePrefix ./work/star/${id}/${id}. \
--outSAMtype None

#第二轮
STAR --runThreadN 12 \
--genomeDir ${ref_idx_star} \
--readFilesIn ./work/fq/${id}/${id}.fastq \
--outFileNamePrefix ./work/star/${id}/${id}. \
--sjdbFileChrStartEnd ./work/star/${id}/*SJ.out.tab \
--sjdbGTFfile ${ref_gtf} \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM XS \
--quantMode TranscriptomeSAM