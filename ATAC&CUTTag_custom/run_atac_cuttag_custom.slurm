#!/bin/bash
# This is a custom script to run ATACseq_pipeline, version 4.
#     ref. https://www.encodeproject.org/atac-seq/
#          https://github.com/kundajelab/atac_dnase_pipelines          
#     dependencies & softwares not independent yet
#         still have chance to meet problems
# Shaorui, Liu ; 2019.02.28 in THU LLab; liushaorui@mail.bnu.edu.cn. 
#
# modified to fit CUTTag data processing on WLSC clusters.
#     ref. https://yezhengstat.github.io/CUTTag_tutorial    
# Shaorui, Liu ; 2021,05,27 in Westlake XuLab; liushaorui@mail.bnu.edu.cn.
#                2021.09 modified for new Slurm system.

#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH -J seq_name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 12
#SBATCH --mem-per-cpu=6GB
#SBATCH -o %j_slurm.out
#SBATCH -e %j_slurm.err

# env variable
# only set conda/bin into $PATH
#   all bioinfo tools/softwares used in this script were installed with conda
source /home/xuhepingLab/liushaorui/.bashrc

# first build a sample list, specify the path
# list_sample is like:
# sample_1
# sample_2
# sample_3
list_sample=''

# specify input & output directory
#     fastq file must be like:  sample_R1.fastq.gz/sample_R2.fastq.gz
#         or have to manually build soft links with matched file names     
fastq_dir=''

#output
#  usually the working directory
output_dir=''


# specify Threads then others below could follow
Threads='12'
Threads_cut=$(echo ${Threads}/2-1 |bc)    # threads for cutadapt each fastq, each get half and would keep one for file processing

java_memory=70G

## specify Org, option available at present : mm10, hg38 , others get error ! and will exit .

Org='mm10'
#Org='hg38'

## then an easy judgement will choose Org related annotaions automatically
#
if  [ $Org == 'hg38'  ]
then
Index='/storage/xuhepingLab/0.share/genomics/human/encode/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta'
blacklist='/storage/xuhepingLab/0.share/genomics/human/encode/hg38.blacklist.bed'
filter1='v'
filter2='V'
TSS_BED='/storage/xuhepingLab/0.share/genomics/human/encode/hg38_gencode_tss_unique.bed'
TSS_extend='/storage/xuhepingLab/0.share/genomics/human/encode/hg38_gencode_tss_unique.ext2k.bed'
CHROMSIZES='/storage/xuhepingLab/0.share/genomics/human/encode/hg38.filt.chrom.sizes'

RSeQC_model='/storage/xuhepingLab/0.share/genomics/human/RSeQC/hg38_Gencode_V28.bed'
RSeQC_hk='/storage/xuhepingLab/0.share/genomics/human/RSeQC/hg38.HouseKeepingGenes.bed'

elif  [ $Org == 'mm10'  ]
then
Index='/storage/xuhepingLab/0.share/genomics/mouse/encode/mm10_no_alt_analysis_set_ENCODE.fasta'
blacklist='/storage/xuhepingLab/0.share/genomics/mouse/encode/mm10.blacklist.bed'
filter1='random'
filter2='chrUn'
TSS_BED='/storage/xuhepingLab/0.share/genomics/mouse/encode/mm10_gencode_tss_unique.bed'
TSS_extend='/storage/xuhepingLab/0.share/genomics/mouse/encode/mm10_gencode_tss_unique.ext2k.bed'
CHROMSIZES='/storage/xuhepingLab/0.share/genomics/mouse/encode/mm10.filt.chrom.sizes'

RSeQC_model='/storage/xuhepingLab/0.share/genomics/mouse/RSeQC/mm10_Gencode_VM18.bed'
RSeQC_hk='/storage/xuhepingLab/0.share/genomics/mouse/RSeQC/mm10.HouseKeepingGenes.bed'

else
echo "Org option wrong, script exit, please check again !"
exit 1
fi
#

tss_plot='/storage/xuhepingLab/0.share/pipelines/common_code/kundaje_atac/tss_plot.py'
picard_jar='/storage/xuhepingLab/0.share/applications/miniconda3/share/picard-2.23.8-0/picard.jar'

# major script

for sample in $(cat ${list_sample})
do

####  pre- process

align_dir=$(echo ${output_dir}/${sample}/align )
qc_dir=$(echo  ${output_dir}/${sample}/qc )
mkdir -p  ${align_dir}
mkdir -p  ${qc_dir}

##  detect  adaptor
# use python script in Kundaje pipeline, which was from https://github.com/nboley/GGR_code
#                                              but link unreachable now

python3 \
/storage/xuhepingLab/0.share/pipelines/common_code/detect_adapter.py \
${fastq_dir}/${sample}_R1.fastq.gz \
> ${qc_dir}/adaptor1.log

python3 \
/storage/xuhepingLab/0.share/pipelines/common_code/detect_adapter.py \
${fastq_dir}/${sample}_R2.fastq.gz \
> ${qc_dir}/adaptor2.log

adaptor_seq1=$(cat ${qc_dir}/adaptor1.log |sed -n 9p |cut -f 3 )
adaptor_seq2=$(cat ${qc_dir}/adaptor2.log |sed -n 9p |cut -f 3 )

## 
# cutadapt

cutadapt \
-j ${Threads_cut} \
-m 20 -e 0.1 -O 3 \
-q 30 --quality-base=33 \
-a ${adaptor_seq1} \
-A ${adaptor_seq2} \
-o ${align_dir}/${sample}_R1.trimmed.fastq.gz \
-p ${align_dir}/${sample}_R2.trimmed.fastq.gz \
${fastq_dir}/${sample}_R1.fastq.gz \
${fastq_dir}/${sample}_R2.fastq.gz \
> ${qc_dir}/${sample}_cutadapt_report.txt

## fastqc

fastqc \
-t ${Threads_cut} \
${fastq_dir}/${sample}_R1.fastq.gz \
${fastq_dir}/${sample}_R2.fastq.gz \
-o ${qc_dir} \
-d ${qc_dir}/

fastqc \
-t ${Threads_cut} \
${align_dir}/${sample}_R1.trimmed.fastq.gz \
${align_dir}/${sample}_R2.trimmed.fastq.gz \
-o ${qc_dir} \
-d ${qc_dir}/


#### alignment
##   memory require : one  thread for bowtie2 is 5G

# about the alignment parameter -t -q -N1 -L 25
# could use --local or --end-to-end to replace

bowtie2 \
-I 10 \
-X 2000 \
--mm \
--local --very-sensitive-local \
--no-mixed --no-discordant \
--threads ${Threads} \
-x ${Index} \
-1 ${align_dir}/${sample}_R1.trimmed.fastq.gz \
-2 ${align_dir}/${sample}_R2.trimmed.fastq.gz \
2>${qc_dir}/${sample}_bowtie2.log |\
samtools view -@ ${Threads} -Su /dev/stdin |\
samtools sort -@ ${Threads} -m 4G - > ${align_dir}/${sample}.bam

samtools index -@ ${Threads} ${align_dir}/${sample}.bam

#### post-align

## filter 1/3    
# (use MAPQ instead of processing uniquely mapped reads;  uniquely mapping rarely mentioned today )
# flag: filter 1804=1024+512+256+8+4 ; get 2
# MAPQ > 30
# sort by name 
samtools view -F 1804 -f 2 -q 30 -@ ${Threads} -u ${align_dir}/${sample}.bam |\
samtools sort -@ ${Threads} -m 4G -n /dev/stdin -o ${align_dir}/${sample}.tmp.filt.bam

## filter 2/3
# fix mate info of name sorted bam
# sort by coordinate again 
samtools fixmate -@ ${Threads} -r ${align_dir}/${sample}.tmp.filt.bam ${align_dir}/${sample}.tmp.fixmate.bam
samtools view -F 1804 -f 2 -@ ${Threads} -u ${align_dir}/${sample}.tmp.fixmate.bam |\
samtools sort -@ ${Threads} -m 4G /dev/stdin -o ${align_dir}/${sample}.filter.bam

# just need filter.bam for next ..
rm ${align_dir}/${sample}.tmp.filt.bam
rm ${align_dir}/${sample}.tmp.fixmate.bam

## filter 3/3
# picard mark duplicates (not remove) 
# use samtools view -F 1024(in 1804) to filter, better than picard ?

# had better specify a java temp path for Markduplicates 
#     or it might cause error when the default system path is full
java_temp=$(echo ${output_dir}/java_temp )
mkdir -p ${java_temp}

java -Xmx${java_memory} -Djava.io.tmpdir=${java_temp} \
-jar ${picard_jar} MarkDuplicates \
INPUT=${align_dir}/${sample}.filter.bam \
OUTPUT=${align_dir}/${sample}.dupmark.bam \
METRICS_FILE=${qc_dir}/${sample}.dup.qc \
VALIDATION_STRINGENCY=LENIENT \
ASSUME_SORTED=true \
REMOVE_DUPLICATES=false

samtools view -F 1804 -f 2 -@ ${Threads} -b -u ${align_dir}/${sample}.dupmark.bam |\
samtools sort -@ ${Threads} -m 4G /dev/stdin -o ${align_dir}/${sample}.nodup.bam

samtools index -@ ${Threads} ${align_dir}/${sample}.dupmark.bam
samtools index -@ ${Threads} ${align_dir}/${sample}.nodup.bam

## add one more step to filter blacklist
#     dupmark here only filter for blacklist, no 'samtools view -F -f', or get dupmark.clean as same as nodup.clean
bedtools intersect -v -abam ${align_dir}/${sample}.dupmark.bam -b ${blacklist} |\
samtools view -@ ${Threads} -S -h -b |\
samtools sort -@ ${Threads} -m 4G  /dev/stdin -o  ${align_dir}/${sample}.dupmark.clean.bam

samtools index -@ ${Threads} ${align_dir}/${sample}.dupmark.clean.bam


bedtools intersect -v -abam ${align_dir}/${sample}.nodup.bam -b ${blacklist} |\
samtools view -F 1804 -f 2 -@ ${Threads} -S -h -b |\
samtools sort -@ ${Threads} -m 4G  /dev/stdin -o  ${align_dir}/${sample}.nodup.clean.bam

samtools index -@ ${Threads} ${align_dir}/${sample}.nodup.clean.bam


# flagstat ?
samtools flagstat -@ ${Threads} ${align_dir}/${sample}.bam > ${qc_dir}/${sample}.flagstat
samtools flagstat -@ ${Threads} ${align_dir}/${sample}.dupmark.bam > ${qc_dir}/${sample}.dupmark.flagstat
samtools flagstat -@ ${Threads} ${align_dir}/${sample}.dupmark.clean.bam > ${qc_dir}/${sample}.dupmark.clean.flagstat
samtools flagstat -@ ${Threads} ${align_dir}/${sample}.nodup.bam > ${qc_dir}/${sample}.nodup.flagstat
samtools flagstat -@ ${Threads} ${align_dir}/${sample}.nodup.clean.bam > ${qc_dir}/${sample}.nodup.clean.flagstat

# library complexity
#    use dupmark.bam
echo "TotalPair,DictinctPair,OnePair,TwoPair,NRF=Distinct/Total,PBC1=OnePair/Distinct,PBC2=OnePair/TwoPair" >  ${qc_dir}/${sample}.pbc_qc.csv
bedtools bamtobed -i ${align_dir}/${sample}.dupmark.bam |\
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' |\
grep -v 'chrM' |sort |uniq -c |\
awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} \
{m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; \
printf "%d,%d,%d,%d,%f,%f,%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' >>  ${qc_dir}/${sample}.pbc_qc.csv

# bigwig
# extend
bamCoverage \
-p ${Threads} \
-e 250 \
-bs 100 \
--normalizeUsing RPKM \
-b ${align_dir}/${sample}.dupmark.clean.bam \
-o ${align_dir}/${sample}.dupmark.rpkm.bw

bamCoverage \
-p ${Threads} \
-e 250 \
-bs 100 \
--normalizeUsing RPKM \
-b ${align_dir}/${sample}.nodup.clean.bam \
-o ${align_dir}/${sample}.nodup.rpkm.bw

# remove temp files
#    or just hold for some time , may check for some possible usage
rm ${align_dir}/${sample}.filter.bam
#rm ${align_dir}/${sample}.dupmark.bam               # CUTTag may not need to remove duplicates

## nodup.bam to insert.bed
# filter with blacklist
# bam to bed ,and filter ,and shift 
# get insert.bed
# get insert.ext.bed , ext50 & ext250, abd
# 
# delete additional bed file processing, 
#   could do them manually if needed


# input clean.bam, so no need to filter chrM and chrUn again
bedtools bamtobed -i ${align_dir}/${sample}.dupmark.clean.bam |\
gawk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' |\
gawk -F "\t" 'BEGIN{OFS=FS}{if($6=="+"){$2=$2+4} else if($6=="-"){$3=$3-5} print $0}' \
> ${align_dir}/${sample}.dupmark.shift.bed

cat ${align_dir}/${sample}.dupmark.shift.bed |gawk '{if($6=="+"){print $1"\t"$2-1"\t"$2;} else{print $1"\t"$3-1"\t"$3}}' \
> ${align_dir}/${sample}.dupmark.insert.bed


bedtools bamtobed -i ${align_dir}/${sample}.nodup.clean.bam |\
gawk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' |\
gawk -F "\t" 'BEGIN{OFS=FS}{if($6=="+"){$2=$2+4} else if($6=="-"){$3=$3-5} print $0}' \
> ${align_dir}/${sample}.nodup.shift.bed

cat ${align_dir}/${sample}.nodup.shift.bed |gawk '{if($6=="+"){print $1"\t"$2-1"\t"$2;} else{print $1"\t"$3-1"\t"$3}}' \
> ${align_dir}/${sample}.nodup.insert.bed


#### add other qc like library complexity
# see above
## fragment_size histogram
# could do better plot for this distribution in R manually if needed

java -Xmx${java_memory} -jar ${picard_jar} CollectInsertSizeMetrics \
I=${align_dir}/${sample}.dupmark.clean.bam \
O=${qc_dir}/${sample}.dupmark.fragsize.txt \
H=${qc_dir}/${sample}.dupmark.fragsize.pdf \
INCLUDE_DUPLICATES=true \
VERBOSITY=ERROR QUIET=TRUE \
W=2000

java -Xmx${java_memory} -jar ${picard_jar} CollectInsertSizeMetrics \
I=${align_dir}/${sample}.nodup.clean.bam \
O=${qc_dir}/${sample}.nodup.fragsize.txt \
H=${qc_dir}/${sample}.nodup.fragsize.pdf \
VERBOSITY=ERROR QUIET=TRUE \
W=2000


## tss
# load conda env built for kundaje pipeline
#     dependency in this python script hasn't been independent yet
# tss_plot.py is a simplified version generated from the powful QC part of their pipeline
# method in this python script is the one used in Jason Nature2013/2015
#     which considers the nearby background signals      
#     (but our current blacklist is only from encode, not the exact one used in the raw paper)
source activate atac_py2

OUTDIR=${qc_dir}/

OUTPREFIX=${sample}.dupmark
FINALBAM=${align_dir}/${sample}.dupmark.clean.bam

python ${tss_plot} \
--outdir $OUTDIR \
--outprefix $OUTPREFIX \
--threads ${Threads} \
--tss ${TSS_BED} \
--finalbam $FINALBAM \
--chromsizes $CHROMSIZES

intersectBed -a ${TSS_extend} -b ${FINALBAM} |wc -l > ${qc_dir}/${sample}.dupmark_reads_in_tss.txt
intersectBed -a ${TSS_extend} -b ${FINALBAM} -wa |sort -u |wc -l > ${qc_dir}/${sample}.dupmark_reads_catched_tss.txt


OUTPREFIX=${sample}.nodup
FINALBAM=${align_dir}/${sample}.nodup.clean.bam

python ${tss_plot} \
--outdir $OUTDIR \
--outprefix $OUTPREFIX \
--threads ${Threads} \
--tss ${TSS_BED} \
--finalbam $FINALBAM \
--chromsizes $CHROMSIZES

intersectBed -a ${TSS_extend} -b ${FINALBAM} |wc -l > ${qc_dir}/${sample}.nodup_reads_in_tss.txt
intersectBed -a ${TSS_extend} -b ${FINALBAM} -wa |sort -u |wc -l > ${qc_dir}/${sample}.nodup_reads_catched_tss.txt

conda deactivate


## summary qc
# note that all reads/fragments here mean one single cutting event
#   there are two kinds of ways to deal with Tn5-based seq data like ATAC
#     one is using real fragments for downstream enrichment
#     the other is only considering the single cutting site
#   we use the second one, so always treat 5' or 3' as independant cutting event then

echo -e "Index,Reads,Reads_mt,Reads_mt%,Reads_noMT,Frag_dup,Frag_dup_mt,Frag_dup_mt%,Frag_dup_noMT,Frag_nodup,Frag_nodup_mt,Frag_nodup_mt%,Frag_nodup_noMT,Alignment_rate,dup%,genomic_dup%,Tss_enrich_dup,TSS_2kb_reads_dup,inTSS_ratio_dup,TSS_2kb_catched_dup,catchTSS_ratio_dup,Tss_enrich_nodup,TSS_2kb_reads_nodup,inTSS_ratio_nodup,TSS_2kb_catched_nodup,catchTSS_ratio_nodup" \
> ${qc_dir}/${sample}.summary.csv

Reads=$(echo $(samtools view -@ ${Threads} ${align_dir}/${sample}.bam |wc -l )/2 |bc)
Reads_m=$(echo $(samtools view -@ ${Threads} ${align_dir}/${sample}.bam |grep chrM |wc -l)/2 |bc)
Reads_m_r=$(echo "scale=4;${Reads_m}/${Reads}" |bc)

Frag_dup=$(echo $(samtools view -@ ${Threads} ${align_dir}/${sample}.dupmark.bam |wc -l)/2 |bc)
Frag_dup_m=$(echo $(samtools view -@ ${Threads} ${align_dir}/${sample}.dupmark.bam |grep chrM |wc -l)/2 |bc)
Frag_dup_m_r=$(echo "scale=4;${Frag_dup_m}/${Frag_dup}" |bc)
Frag_dup_n=$(echo "${Frag_dup}-${Frag_dup_m}" |bc)

Frag_nodup=$(echo $(samtools view -@ ${Threads} ${align_dir}/${sample}.nodup.bam |wc -l)/2 |bc)
Frag_nodup_m=$(echo $(samtools view -@ ${Threads} ${align_dir}/${sample}.nodup.bam |grep chrM |wc -l)/2 |bc)
Frag_nodup_m_r=$(echo "scale=4;${Frag_nodup_m}/${Frag_nodup}" |bc)
Frag_nodup_n=$(echo "${Frag_nodup}-${Frag_nodup_m}" |bc)

Align=$(cat ${qc_dir}/${sample}_bowtie2.log |grep "alignment rate" |gawk '{print $1}' )
Align_d=$(echo 0.01*${Align} |cut -d "%" -f 1 |bc) 

Reads_n=$(printf "%.0f\n" `echo "${Reads}*${Align_d}-${Reads_m}" |bc`)
dup_r=$(echo "scale=4;1-${Frag_nodup}/(${Frag_dup})" |bc )
genomic_dup_r=$(echo "scale=4;(1-${Frag_nodup_n}/${Frag_dup_n})" |bc )

TSS_n=$(cat ${TSS_extend} |wc -l)

TSS_dup=$(printf $(cat ${qc_dir}/${sample}.dupmark_tss-enrich.txt )"\n" )
inTSS_dup=$(cat ${qc_dir}/${sample}.dupmark_reads_in_tss.txt )
inTSS_dup_r=$(echo "scale=4;1/2*${inTSS_dup}/${Frag_dup_n}" |bc)
catchTSS_dup=$(cat ${qc_dir}/${sample}.dupmark_reads_catched_tss.txt )
catchTSS_dup_r=$(echo "scale=4;${catchTSS_dup}/${TSS_n}" |bc)

TSS_nodup=$(printf $(cat ${qc_dir}/${sample}.nodup_tss-enrich.txt )"\n" )
inTSS_nodup=$(cat ${qc_dir}/${sample}.nodup_reads_in_tss.txt )
inTSS_nodup_r=$(echo "scale=4;1/2*${inTSS_nodup}/${Frag_nodup_n}" |bc)
catchTSS_nodup=$(cat ${qc_dir}/${sample}.nodup_reads_catched_tss.txt )
catchTSS_nodup_r=$(echo "scale=4;${catchTSS_nodup}/${TSS_n}" |bc)

echo -e ${sample}","${Reads}","${Reads_m}","${Reads_m_r}","${Reads_n}","${Frag_dup}","${Frag_dup_m}","${Frag_dup_m_r}","${Frag_dup_n}","${Frag_nodup}","${Frag_nodup_m}","${Frag_nodup_m_r}","${Frag_nodup_n}","${Align}","${dup_r}","${genomic_dup_r}","${TSS_dup}","${inTSS_dup}","${inTSS_dup_r}","${catchTSS_dup}","${catchTSS_dup_r}","${TSS_nodup}","${inTSS_nodup}","${inTSS_nodup_r}","${catchTSS_nodup}","${catchTSS_nodup_r} \
>> ${qc_dir}/${sample}.summary.csv

done

#### summary
##

cd ${output_dir}

mkdir bam_file
mkdir bigwig_file
mkdir frag_size
mkdir tss_enrich

mkdir multiqc_all
mkdir fastqc_trim


for sample in $(cat ${list_sample})
do

mv ./${sample}/qc/*R1_fastqc.zip ./multiqc_all/
mv ./${sample}/qc/*R2_fastqc.zip ./multiqc_all/
cp ./${sample}/qc/*cutadapt_report.txt ./multiqc_all/
cp ./${sample}/qc/*bowtie2.log ./multiqc_all/

mv ./${sample}/qc/*R1.trimmed_fastqc.zip ./fastqc_trim/
mv ./${sample}/qc/*R2.trimmed_fastqc.zip ./fastqc_trim/

python /storage/xuhepingLab/0.share/applications/miniconda3/bin/read_distribution.py \
-r ${RSeQC_model} \
-i ./${sample}/align/${sample}.dupmark.clean.bam \
> ./multiqc_all/${sample}.dup.read_distribution.txt

python /storage/xuhepingLab/0.share/applications/miniconda3/bin/read_distribution.py \
-r ${RSeQC_model} \
-i ./${sample}/align/${sample}.nodup.clean.bam \
> ./fastqc_trim/${sample}.nodup.read_distribution.txt

done

cd multiqc_all
multiqc ./ -n multiqc.report

cd ../fastqc_trim
multiqc ./ -n multiqc_trim.report


#
cd ${output_dir}


## summary
echo -e "Index,Reads,Reads_mt,Reads_mt%,Reads_noMT,Frag.dup,Frag.dup_mt,Frag.dup_mt%,Frag.dup_noMT,Frag.nodup,Frag.nodup_mt,Frag.nodup_mt%,Frag.nodup_noMT,Alignment_rate,dup%,genomic_dup%,Tss_enrich.dup,TSS_2kb_reads.dup,inTSS_ratio.dup,TSS_2kb_catched.dup,catchTSS_ratio.dup,Tss_enrich.nodup,TSS_2kb_reads.nodup,inTSS_ratio.nodup,TSS_2kb_catched.nodup,catchTSS_ratio.nodup" \
> summary.csv

for sample in $(cat ${list_sample})
do

#ln -s ${output_new_adv}/${sample}/align/${sample}.*.clean.bam* ./bam_file/
#ln -s ${output_new_adv}/${sample}/align/${sample}.*.rpkm.bw ./bigwig_file/
#ln -s ${output_new_adv}/${sample}/qc/${sample}.*.fragsize.pdf ./frag_size/

ln -s ../${sample}/align/${sample}.*.clean.bam* ./bam_file/
mv ./${sample}/align/${sample}.*.rpkm.bw ./bigwig_file/
mv ./${sample}/qc/${sample}.*.fragsize.pdf ./frag_size/
mv ./${sample}/qc/*.png ./tss_enrich/

cat ${output_dir}/${sample}/qc/${sample}.summary.csv |tail -n 1 >> summary.csv
done



#



