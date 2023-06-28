## Version: vM25  20200429 update
# http://asia.ensembl.org/Mus_musculus/Info/Index
genome_dir='folder-to-store-genome-files'
cd ${genome_dir}

## GRCm38_vM25
# genome.fa 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz

# gtf
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz

# lnRNA gtf
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.long_noncoding_RNAs.gtf.gz

# rRNA ref file from broad
# https://software.broadinstitute.org/cancer/cga/rnaseqc_download
wget https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/rnaseqc/mouse_rRNA.tar.gz

# unzip
# bwa count rRNA counts, need index, and different version not match
# rebuild with bwa-0.7.17
rm mouse_all_rRNA.fasta.*
bwa index mouse_all_rRNA.fasta

# isoform map  :  gene_name"\t"transcript_id
zcat gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz |\
grep -P '\ttranscript\t' |\
awk -F '[\t|\"|;]' '{ print $19"\t"$13 }' \
> gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gene_iso_map
 
## STAR index
# specify folder var
STAR_index=${genome_dir}/STAR_index

# code to submit
STAR \
--runMode genomeGenerate \
--runThreadN 5 \
--genomeDir ${STAR_index} \
--genomeFastaFiles ${genome_dir}/GRCm38.p6.genome.fa \
--sjdbGTFfile ${genome_dir}/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf \
--sjdbOverhang 150

## RSEM index
RSEM_index=${genome_dir}/STAR_index
cd ${RSEM_index}

samtools faidx ${genome_dir}/GRCm38.p6.genome.fa

mkdir GRCm38.p6.genome
ln -s ${genome_dir}/GRCm38.p6.genome.fa ./GRCm38.p6.genome/
ln -s ${genome_dir}/GRCm38.p6.genome.fa.fai ./GRCm38.p6.genome/

# code to submit
rsem-prepare-reference  \
--gtf ${genome_dir}/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf \
--transcript-to-gene-map ${genome_dir}/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gene_iso_map \
-p 5 \
${RSEM_index}/GRCm38.p6.genome \
${RSEM_index}/rsem_trans_index

## gtf details
gtf_detail=${genome_dir}/gtf_detail
cd ${gtf_detail}

cat ${genome_dir}/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf |\
grep -P '\ttranscript\t' |\
awk -F '[\t|\"|;]' '{print $16"\t"$19}' |\
sort -k2 -u > trans_genetype.txt

cat trans_genetype.txt | cut -f 1 |sort |uniq -c |sort -k1,1nr > genecode.vM25.transcript.summary

# protein coding gene list
cat trans_genetype.txt |grep 'protein_coding' |cut -f 2 > list_pc

# only level 1,2
#   check ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/_README.TXT
cat ${genome_dir}/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf |\
awk '{if($0~"level (1|2);"){print $0}}' |\
grep -P '\ttranscript\t' |\
awk -F '[\t|\"|;]' '{print $16"\t"$19}' |\
sort -k2 -u > trans_genetype.lv1_2.txt

cat trans_genetype.lv1_2.txt | cut -f 1 |sort |uniq -c |sort -k1,1nr > genecode.vM25.transcript.lv1_2.summary

cat trans_genetype.lv1_2.txt |grep 'protein_coding' |cut -f 2 > list_pc.lv1_2

## RSeQC scripts and bed files 
# https://rseqc.sourceforge.net/#download

## detect adapter
# https://github.com/kundajelab/atac_dnase_pipelines/raw/master/utils/detect_adapter.py

#### end ####

