##########################################################################
# RNA-seq Time Course: PPD-b stimulated vs unstimulated peripheral blood #
# paired-end reads.                                                      #
#       --- Linux bioinformatics workflow for known sense genes ---      #
##########################################################################

# Based on the pipeline created by Nalpas, N.C. (2014) 
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474
# Author of current version (4.0.0): Correia, C.N.
# DOI badge of current version:
# Last updated on: 30/07/2017

################################
# Download and files check sum #
################################

# All files were downloaded from MSU in 2013. At the time, they were
# also md5sum checked and renamed.

# File names correspond to:
# AnimalNumber_TimePoint_AnimalGroup_PairedEndTag_LaneNumber_fastq.gz

###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering \
--noextract --nogroup -t 2 \
/workspace/storage/kmcloughlin/RNAseqTimeCourse/A6511_W10_P_R1_001.fastq.gz

### Moved this folder to my laptop via SCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering $file" \
>> fastqc.sh; done;

# Split and run all scripts on Stampede:
split -d -l 70 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all the files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp; \
done;

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; done

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -rf $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp

##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html

# Create a working directory for filtered reads:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/
cd $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/

# Run ngsShoRT in one pair of reads to check if it's working:
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 /workspace/storage/kmcloughlin/RNAseqTimeCourse/A6522_W10_U_R1_004.fastq.gz \
-pe2 /workspace/storage/kmcloughlin/RNAseqTimeCourse/A6522_W10_U_R2_004.fastq.gz \
-o /home/ccorreia/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6522_W10_U_004 \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

# Create bash scripts to perform filtering of each FASTQ file, keeping the
# sequencing lane information:
for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/001/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.001.sh; done;

for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_002.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/002/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.002.sh; done;

for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_003.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/003/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.003.sh; done;

for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_004.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/004/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.004.sh; done;

for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_005.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/005/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.005.sh; done;

# Run all scripts on Stampede:
for script in `ls filtering.00*.sh`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check that all files were processed:
for file in `ls filtering.00*.sh.nohup`; \
do grep -o 'Done-MAIN' $file | wc -l; done

# Compress files with discarded reads:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence \
-name extracted_*.txt`; do echo "gzip -9 $file" >> discarded_compression.sh; \
done;

# Split and run all scripts on Stampede:
split -d -l 70 discarded_compression.sh discarded_compression.sh.
for script in `ls discarded_compression.sh.*`
do
chmod 755 $script
nohup ./$script &
done

################################################
# FastQC quality check of filtered FASTQ files #
################################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering
cd $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering \
--noextract --nogroup -t 10 \
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6644_W6_F_001/trimmed_A6644_W6_F_R1_001.fastq.gz

### Moved this folder to my laptop via SCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all filtered
# FASTQ files:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/ \
-name *_R*_00*.fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering $file" \
>> fastqc_filt.sh; done;

# Split and run all scripts on Stampede:
split -d -l 75 fastqc_filt.sh fastqc_filt.sh.
for script in `ls fastqc_filt.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all the files were processed:
for file in `ls fastqc_filt.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp; \
done;

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp \
-name summary.txt`; do more $file >> reports_post-filtering.txt; done

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post-filtering.txt; \
done

# Remove temporary folder and its files:
rm -rf $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp

##############################################################################
# Alignment of FASTQ files against the Bos taurus reference genome with STAR #
##############################################################################

# Required software is STAR 2.5.1b, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# NCBI changed the location of several genome and annotation files in their FTP
# server on December 2016, which means that the original links listed below
# do not work anymore. The same files used in this analysis can now be found at:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1

# Download Bos taurus reference genome, version UMD3.1.1 from NCBI:
mkdir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/source_file
cd /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/source_file
nohup wget nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz &

gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz

# Download annotation file for UMD3.1.1 NCBI Bos taurus genomic annotation (GCF_000003055.6):
mkdir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file
cd /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz
gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz

# Generate genome indexes files using annotations:
mkdir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index
cd /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index

nohup STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index \
--genomeFastaFiles \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/source_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna \
--sjdbGTFfile /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 \
--outFileNamePrefix \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index/Btau-UMD3.1.1 &

# Create and enter alignment working directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment
cd $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment

# Mapping reads from one FASTQ file to the indexed genome,
# to check if it works well:
nohup STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index/ \
--readFilesIn \
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_001/trimmed_A6511_W10_F_R1_001.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_002/trimmed_A6511_W10_F_R1_002.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_003/trimmed_A6511_W10_F_R1_003.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_004/trimmed_A6511_W10_F_R1_004.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_005/trimmed_A6511_W10_F_R1_005.fastq.gz \
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_001/trimmed_A6511_W10_F_R2_001.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_002/trimmed_A6511_W10_F_R2_002.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_003/trimmed_A6511_W10_F_R2_003.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_004/trimmed_A6511_W10_F_R2_004.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_005/trimmed_A6511_W10_F_R2_005.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./A6511_W10_F_ \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &

# Create a bash script to perform alignment of paired FASTQ files:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence \
-name *_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(_00.)/_002/g'`; \
file3=`echo $file | perl -p -e 's/(_00.)/_003/g'`; \
file4=`echo $file | perl -p -e 's/(_00.)/_004/g'`; \
file5=`echo $file | perl -p -e 's/(_00.)/_005/g'`; \
read1=`echo $file | perl -p -e 's/(R1_00.)/R2_001/'`; \
read2=`echo $file2 | perl -p -e 's/(R1_00.)/R2_002/'`; \
read3=`echo $file3 | perl -p -e 's/(R1_00.)/R2_003/'`; \
read4=`echo $file4 | perl -p -e 's/(R1_00.)/R2_004/'`; \
read5=`echo $file5 | perl -p -e 's/(R1_00.)/R2_005/'`; \
sample=`basename $file | perl -p -e 's/\_R1_001\.fastq\.gz//'`; \
foldername=`basename $sample | perl -p -e 's/trimmed\_//'`; \
echo "mkdir $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment/$foldername; \
cd $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment/$foldername; \
STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index/ \
--readFilesIn $file,$file2,$file3,$file4,$file5 \
$read1,$read2,$read3,$read4,$read5 --readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./${foldername}_ --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 70 alignment.sh alignment.sh.
for script in `ls alignment.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check nohup.out file to see how many jobs finished successfully:
grep -c 'Finished successfully' alignment.sh.00.nohup
grep -c 'Finished successfully' alignment.sh.01.nohup

# Merge all STAR log.final.out files into a single file:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment \
-name *Log.final.out`; \
do perl /home/nnalpas/SVN/star_report_opener.pl -report $file; done;

#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment
cd $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment

# Create a bash script to perform FastQC quality check on aligned SAM files:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment $file" >> \
fastqc_aligned.sh; done;

# Split and run all scripts on Stampede
split -d -l 35 fastqc_aligned.sh fastqc_aligned.sh.
for script in `ls fastqc_aligned.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Delete all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment/tmp; \
done

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment/tmp \
-name summary.txt`; do more $file >> reports_post-alignment.txt; done

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Check if all files were processed:
grep -c '##FastQC' basic_stats_post_alignment.txt
grep -c 'Basic Statistics' reports_post-alignment.txt
grep -c 'Analysis complete' fastqc_aligned.sh.00.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.01.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.02.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.03.nohup

# Remove temporary folder:
rm -r tmp/

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################

# Required package is featureCounts, which is part of Subread 1.5.1 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directories:
cd $HOME/scratch/PPDbRNAseqTimeCourse/
mkdir -p Count_summarisation/sense
cd $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense

# Run featureCounts with one sample to check if it is working fine:
featureCounts -a \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
-B -p -C -R -s 1 -T 15 -t gene -g Dbxref -o ./counts.txt \
$HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment/A6511_W10_F/A6511_W10_F_Aligned.out.bam

# Create a bash script to run featureCounts on BAM file containing multihits and
# uniquely mapped reads using the stranded parameter:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment \
-name *_Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Aligned.out.bam//'`; \
echo "mkdir $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/$sample; \
cd $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/$sample; \
featureCounts -a \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
-B -p -C -R -s 1 -T 10 -t gene -g Dbxref \
-o ${sample}_sense-counts.txt $file" >> sense_count.sh; done

# Split and run all scripts on Stampede:
split -d -l 70 sense_count.sh sense_count.sh.
for script in `ls sense_count.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
grep -c 'Read assignment finished.' sense_count.sh.00.nohup
grep -c 'Read assignment finished.' sense_count.sh.01.nohup


# Create bash script to merge stats info from .featureCounts from all samples
# into a single file:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/ \
-name *.featureCounts`; do echo echo \
"\`basename $file\` \`cut $file -f2 | sort | uniq -c | perl -p -e 's/\n/ /'\` >> \
annotation_summary_sense.txt" >> annotation_summary_sense.sh
done

# Split and run scripts on Stampede:
split -d -l 70 annotation_summary_sense.sh annotation_summary_sense.sh.
for script in `ls annotation_summary_sense.sh.*`
do
chmod 755 $script
nohup ./$script &
done

# Check that all files were processed:
grep -c '.featureCounts' annotation_summary_sense.txt

# Copy all *sense-counts.txt files to temporary folder:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/tmp

for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/ \
-name *sense-counts.txt`; do cp $file \
-t $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/tmp; \
done

# Transfer all files from tmp to laptop, then remove tmp folder:
rm -r tmp


########################################
# R analysis of gene counts with edgeR #
########################################

# Subsequent sense genes analyses were performed using the R statistical
# software and the edgeR package.
# Please go to file: 01-PPDb-RNA-seq_paired_sense.R














