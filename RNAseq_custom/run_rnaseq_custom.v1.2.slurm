#!/bin/bash

##
# This is a custom script to run RNAseq_pipeline
#   could run any normal RNAseq without spike-in
#             bulk or single like smartseq2    
#   ref: check https://www.encodeproject.org/pages/pipelines/#RNA-seq
# Shaorui, Liu; 2019; liushaorui@mail.bnu.edu.cn.
#   modified on WLSC, 202009, available for xuhp group
#   modified on WLSC Slurm, 202109, available for Xu&He Group
##  


#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH -J SS2_name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 8
#SBATCH --mem-per-cpu=8GB
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

# matrix name, output will be like RNAseq.${mat_name}.counts.matrix
mat_name='date_SS2_name'

## specify Org, option available at present : mm10, hg38 , others get error ! and will exit .
##     un-annotate one option and annotate the other 
## current version: mm10_vM25 and hg38_v34, Gencode 2020.04.29 updated. 

Org='mm10'
#Org='hg38'

#### general user should stop here ####

## then an easy judgement will choose Org related annotaions automatically
#
if  [ $Org == 'hg38'  ]
then
RRNA_FA='/storage/xuhepingLab/0.share/genomics/human/GRCh38_v34/human_rRNA/human_all_rRNA.fasta'
STAR_INDEX='/storage/xuhepingLab/0.share/genomics/human/GRCh38_v34/STAR_index'
RSEM_TRANS='/storage/xuhepingLab/0.share/genomics/human/GRCh38_v34/RSEM_index/rsem_trans_index'
RSeQC_model='/storage/xuhepingLab/0.share/genomics/human/RSeQC/hg38_Gencode_V28.bed'
RSeQC_hk='/storage/xuhepingLab/0.share/genomics/human/RSeQC/hg38.HouseKeepingGenes.bed'

elif  [ $Org == 'mm10'  ]
then
RRNA_FA='/storage/xuhepingLab/0.share/genomics/mouse/GRCm38_vM25/mouse_rRNA/mouse_all_rRNA.fasta'
STAR_INDEX='/storage/xuhepingLab/0.share/genomics/mouse/GRCm38_vM25/STAR_index'
RSEM_TRANS='/storage/xuhepingLab/0.share/genomics/mouse/GRCm38_vM25/RSEM_index/rsem_trans_index'
RSeQC_model='/storage/xuhepingLab/0.share/genomics/mouse/RSeQC/mm10_Gencode_VM18.bed'
RSeQC_hk='/storage/xuhepingLab/0.share/genomics/mouse/RSeQC/mm10.HouseKeepingGenes.bed'

else
echo "Org option wrong, script exit, please check again !"
exit 1
fi
#

# major script

for sample in $(cat ${list_sample})
do

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
-j 3 -m 20 -e 0.1 -O 3 \
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
-t 6 \
${fastq_dir}/${sample}_R1.fastq.gz \
${fastq_dir}/${sample}_R2.fastq.gz \
-o ${qc_dir} \
-d ${qc_dir}/

fastqc \
-t 6 \
${align_dir}/${sample}_R1.trimmed.fastq.gz \
${align_dir}/${sample}_R2.trimmed.fastq.gz \
-o ${qc_dir} \
-d ${qc_dir}/

## RRNA count 
# used to use custom srcipt
# rewrite it with bwa
# ref: RNA-SeQC from broad 
#   https://github.com/francois-a/rnaseqc  CountAligned.java

mkdir -p ${qc_dir}/RRNA

bwa aln \
-t 8 \
${RRNA_FA} \
${align_dir}/${sample}_R1.trimmed.fastq.gz > ${qc_dir}/RRNA/rRNA_sa1.sai

bwa aln \
-t 8 \
${RRNA_FA} \
${align_dir}/${sample}_R2.trimmed.fastq.gz > ${qc_dir}/RRNA/rRNA_sa2.sai

bwa sampe \
${RRNA_FA} \
${qc_dir}/RRNA/rRNA_sa1.sai \
${qc_dir}/RRNA/rRNA_sa2.sai \
${align_dir}/${sample}_R1.trimmed.fastq.gz \
${align_dir}/${sample}_R2.trimmed.fastq.gz \
> ${qc_dir}/RRNA/rRNA_aln-pe.sam

samtools view -@ 6 -q 30 ${qc_dir}/RRNA/rRNA_aln-pe.sam |\
wc -l > ${qc_dir}/RRNA//rRNA.sam.frag_count

rm ${qc_dir}/RRNA/rRNA_sa1.sai ${qc_dir}/RRNA/rRNA_sa2.sai ${qc_dir}/RRNA/rRNA_aln-pe.sam


#### alignment
##  only star-rsem available at present 

# star
mkdir -p ${align_dir}/STAR
STAR \
--genomeDir ${STAR_INDEX} \
--runThreadN 8 \
--readFilesIn \
${align_dir}/${sample}_R1.trimmed.fastq.gz \
${align_dir}/${sample}_R2.trimmed.fastq.gz \
--readFilesCommand gunzip -c \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--limitBAMsortRAM 50000000000 \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFileNamePrefix ${align_dir}/STAR/ \
--outSAMstrandField intronMotif \
--quantMode TranscriptomeSAM

# bam file process
#     just index & bigwig convertion for igv visualization
samtools index -@ 6  ${align_dir}/STAR/Aligned.sortedByCoord.out.bam

# v1.2
# add CPM as normalization, or the signal would vary from library size
bamCoverage \
--normalizeUsing CPM \
-p 6 \
-b  ${align_dir}/STAR/Aligned.sortedByCoord.out.bam \
-o  ${align_dir}/STAR/Aligned.sortedByCoord.out.bw

# rsem
mkdir -p ${align_dir}/RSEM_STAR
rsem-calculate-expression \
-p 8 \
--no-bam-output \
--paired-end \
--alignments ${align_dir}/STAR/Aligned.toTranscriptome.out.bam \
${RSEM_TRANS} \
${align_dir}/RSEM_STAR/RSEM_STAR_${sample}

done


#### summary
##

cd ${output_dir}

mkdir bigwig_file
mkdir bam_file
mkdir gene_result
mkdir isoform_result
mkdir multiqc_all
mkdir fastqc_trim

for line in $(cat ${list_sample})
do
# put all bigwig (for IGV ), gene result ( just keep col 1,5,6,7 ) into a new dir
mv ./${line}/align/STAR/Aligned.sortedByCoord.out.bw ./bigwig_file/${line}.bw
cat ./${line}/align/RSEM_STAR/*.genes.results |gawk 'OFS="\t"{print $1,$5,$6,$7 }' > ./gene_result/${line}.genes.results
cat ./${line}/align/RSEM_STAR/*.isoforms.results |gawk 'OFS="\t"{print $1,$5,$6,$7 }' > ./isoform_result/${line}.isoforms.results

# multiqc
mv ./${line}/qc/*R1_fastqc.zip ./multiqc_all/
mv ./${line}/qc/*R2_fastqc.zip ./multiqc_all/
cp ./${line}/qc/*cutadapt_report.txt ./multiqc_all/
cp ./${line}/align/STAR/*Log.final.out ./multiqc_all/${line}.Log.final.out
cp ./${line}/align/STAR/*Log.out ./multiqc_all/${line}.Log.out
cp ./${line}/align/STAR/*Log.progress.out ./multiqc_all/${line}.Log.progress.out
cp ./${line}/align/RSEM_STAR/*${line}.stat/*cnt ./multiqc_all/${line}.cnt

mv ./${line}/qc/*R1.trimmed_fastqc.zip ./fastqc_trim/
mv ./${line}/qc/*R2.trimmed_fastqc.zip ./fastqc_trim/
done

# create link to bamfile of each sample (sort by genome coordinate)
for line in $(cat ${list_sample})
do
ln -s ./${line}/align/STAR/Aligned.sortedByCoord.out.bam bam_file/${line}.bam
ln -s ./${line}/align/STAR/Aligned.sortedByCoord.out.bam.bai bam_file/${line}.bam.bai
done

   
## summary.csv
# sample, total reads, uniqulely mapping reads/ratio, multiple mapping reads/ratio, rRNA, rRNA_r= rRNA/(unique+multiple)
# add tpm>0, tpm>2 gene count
# add two types of mapping, 
#    type I:  genome + transcriptome
#    type II:  transcriptome


echo -e "Sample,TotalReads(raw),TotalReads(trimmed),TrimmedRatio,star(geno+trans),Unique,UniqueRatio,Mutiple,MutipleRatio,star(trans),Unique,UniqueRatio,Mutiple,MutipleRatio,rRNA,rRNARatio,genes(TPM>0),genes(TPM>2)" > summary.csv
for line in $(cat ${list_sample})
do
Reads_raw=$(cat ./${line}/qc/*cutadapt_report.txt |grep "Total read pairs" |gawk -F ":" '{print $2}' |sed 's/,//g' |sed 's/ //g')
Reads=$(cat ./${line}/align/STAR/Log.final.out |grep "Number of input reads" |cut -f 2 )
Trim_r=$(echo "scale=2;${Reads}*100/${Reads_raw}" |bc )"%"

Unique=$(cat ./${line}/align/STAR/Log.final.out |grep "Uniquely mapped reads number" |cut -f 2 )
Unique_r=$(cat ./${line}/align/STAR/Log.final.out |grep "Uniquely mapped reads %" |cut -f 2 )
Multiple=$(cat ./${line}/align/STAR/Log.final.out |grep "Number of reads mapped to multiple loci" |cut -f 2 )
Multipe_r=$(cat ./${line}/align/STAR/Log.final.out |grep "% of reads mapped to multiple loci" |cut -f 2 )

Uniq_trans=$(cat ./${line}/align/RSEM_STAR/*${line}.stat/*cnt |sed -n '2p' |cut -d " " -f 1 )
Uniq_trans_r=$(echo "scale=2;${Uniq_trans}*100/${Reads}" |bc )"%"
Multi_trans=$(cat ./${line}/align/RSEM_STAR/*${line}.stat/*cnt |sed -n '2p' |cut -d " " -f 2 )
Multi_trans_r=$(echo "scale=2;${Multi_trans}*100/${Reads}" |bc )"%"

rRNA=$(cat ./${line}/qc/RRNA/rRNA.sam.frag_count )
rRNA_r=$(echo "scale=2;${rRNA}*100/(${Unique}+${Multiple})" |bc )"%"

genes_1=$(cat  gene_result/${line}.genes.results |cut -f 3 |tail -n +2 |gawk '$1>0' |wc -l )
genes_2=$(cat  gene_result/${line}.genes.results |cut -f 3 |tail -n +2 |gawk '$1>2' |wc -l )

echo ${line},${Reads_raw},${Reads},${Trim_r},"",${Unique},${Unique_r},${Multiple},${Multipe_r},"",${Uniq_trans},${Uniq_trans_r},${Multi_trans},${Multi_trans_r},${rRNA},${rRNA_r},${genes_1},${genes_2} >> summary.csv
done


## Matrix
# counts
idx=$(cat ${list_sample} |head -n 1)

for fmt in {"gene","isoform"}
do
cat ${fmt}_result/${idx}.${fmt}s.results | cut -f 1 > RNAseq.${mat_name}.counts.${fmt}.matrix
for line in $(cat ${list_sample})
do
echo ${line} > ${line}.tmp
cat  ${fmt}_result/${line}.${fmt}s.results |cut -f 2 |tail -n +2 >> ${line}.tmp
paste -d "\t" RNAseq.${mat_name}.counts.${fmt}.matrix ${line}.tmp > RNAseq.${mat_name}.counts.${fmt}.matrix.tmp
rm RNAseq.${mat_name}.counts.${fmt}.matrix
mv RNAseq.${mat_name}.counts.${fmt}.matrix.tmp RNAseq.${mat_name}.counts.${fmt}.matrix
rm ${line}.tmp
done
done

# TPM
idx=$(cat ${list_sample} |head -n 1)

for fmt in {"gene","isoform"}
do
cat ${fmt}_result/${idx}.${fmt}s.results | cut -f 1 > RNAseq.${mat_name}.tpm.${fmt}.matrix
for line in $(cat ${list_sample})
do
echo ${line} > ${line}.tmp
cat  ${fmt}_result/${line}.${fmt}s.results |cut -f 3 |tail -n +2 >> ${line}.tmp
paste -d "\t" RNAseq.${mat_name}.tpm.${fmt}.matrix ${line}.tmp > RNAseq.${mat_name}.tpm.${fmt}.matrix.tmp
rm RNAseq.${mat_name}.tpm.${fmt}.matrix
mv RNAseq.${mat_name}.tpm.${fmt}.matrix.tmp RNAseq.${mat_name}.tpm.${fmt}.matrix
rm ${line}.tmp
done
done

# FPKM
idx=$(cat ${list_sample} |head -n 1)

for fmt in {"gene","isoform"}
do
cat ${fmt}_result/${idx}.${fmt}s.results | cut -f 1 > RNAseq.${mat_name}.fpkm.${fmt}.matrix
for line in $(cat ${list_sample})
do
echo ${line} > ${line}.tmp
cat  ${fmt}_result/${line}.${fmt}s.results |cut -f 4 |tail -n +2 >> ${line}.tmp
paste -d "\t" RNAseq.${mat_name}.fpkm.${fmt}.matrix ${line}.tmp > RNAseq.${mat_name}.fpkm.${fmt}.matrix.tmp
rm RNAseq.${mat_name}.fpkm.${fmt}.matrix
mv RNAseq.${mat_name}.fpkm.${fmt}.matrix.tmp RNAseq.${mat_name}.fpkm.${fmt}.matrix
rm ${line}.tmp
done
done

## multiqc
#
#cd multiqc_all
#multiqc ./ -n ${mat_name}.multiqc
#cd ..

#
#cd fastqc_trim
#multiqc ./ -n ${mat_name}.fastqc_trim
#cd ..


##
mkdir -p analysis

#
mv bigwig_file analysis/
mv *matrix analysis/
mv summary.csv analysis/

#
mkdir -p ${output_dir}/qc_new

for line in $(cat ${list_sample})
do
python /storage/xuhepingLab/0.share/applications/miniconda3/bin/read_distribution.py \
-r ${RSeQC_model} \
-i ${output_dir}/bam_file/${line}.bam \
> ${output_dir}/qc_new/${line}.read_distribution.txt

python /storage/xuhepingLab/0.share/applications/miniconda3/bin/geneBody_coverage.py \
-r ${RSeQC_hk} \
-i ${output_dir}/bam_file/${line}.bam \
-o ${output_dir}/qc_new/${line}

done

cd multiqc_all
mv ../qc_new/* ./
multiqc ./


#### end




