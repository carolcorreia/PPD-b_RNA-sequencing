##########################################################################
# RNA-seq Time Course: PPD-b stimulated vs unstimulated peripheral blood #
# paired-end reads.                                                      #
#       --- Linux bioinformatics workflow for known sense genes ---      #
##########################################################################

# Based on the pipeline created by Nalpas, N.C. (2014) 
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474

# Author: Carolina N. Correia
# Last updated on: 17/04/2021

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

# Required software is FastQC v0.11.9, consult manual/tutorial
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
-name *_P_*.fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering $file" \
>> fastqc.sh; done;

for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *_U_*.fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering $file" \
>> fastqc.sh; done;

# Split and run all scripts on Stampede:
split -d -l 80 fastqc.sh fastqc.sh.
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

####################################
# MultiQC for consolidated reports #
####################################

# Required software is MultiQc version 1.9, consult manual/tutorial
# for details: https://multiqc.info/docs/#-20
cd /home/ccorreia/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering
multiqc .

##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html

# Create a working directory for filtered reads:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/
cd !$

# Run ngsShoRT in one pair of reads to check if it's working:
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 /workspace/storage/kmcloughlin/RNAseqTimeCourse/A6522_W10_U_R1_004.fastq.gz \
-pe2 /workspace/storage/kmcloughlin/RNAseqTimeCourse/A6522_W10_U_R2_004.fastq.gz \
-o /home/ccorreia/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6522_W10_U_004 \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

# Create bash scripts to perform filtering of each FASTQ file, keeping the
# sequencing lane information:
A6514_W10_P_R1_005.fastq.gz
for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *P_R1_001.fastq.gz`; \
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

# Required software is FastQC v0.11.9, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering
cd !$

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

####################################
# MultiQC for consolidated reports #
####################################

# Required software is MultiQc version 1.9, consult manual/tutorial
# for details: https://multiqc.info/docs/#-20
cd /home/ccorreia/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering
multiqc .

##############################################################################
# Alignment of FASTQ files against the Bos taurus reference genome with STAR #
##############################################################################

# Required software is STAR 2.7.8a, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Create and enter the working directory:
mkdir /home/ccorreia/scratch/PPDbRNAseqTimeCourse/ARS-UCD1.2
cd !$

# Download and unzip the Bos taurus reference genome ARS-UCD1.2 from NCBI:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9913/106/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz
gunzip GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz

# Download and unzip annotation file for ARS-UCD1.2 NCBI Bos taurus genomic Annotation Release 106:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9913/106/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz
gunzip GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz

# Generate genome indexes files using annotations:
mkdir STAR2.7.8a_index
cd !$

nohup STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir /home/ccorreia/scratch/PPDbRNAseqTimeCourse/ARS-UCD1.2/STAR2.7.8a_index \
--genomeFastaFiles \
/home/ccorreia/scratch/PPDbRNAseqTimeCourse/ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna \
--sjdbGTFfile /home/ccorreia/scratch/PPDbRNAseqTimeCourse/ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 \
--outFileNamePrefix \
/home/ccorreia/scratch/PPDbRNAseqTimeCourse/ARS-UCD1.2/STAR2.7.8a_index/ARS-UCD1.2 &

# Create and enter alignment working directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/STAR2.7.8a_alignment
cd !$

# Mapping reads from one FASTQ file to the indexed genome,
# to check if it works well:
nohup STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir /home/ccorreia/scratch/PPDbRNAseqTimeCourse/ARS-UCD1.2/STAR2.7.8a_index \
--readFilesIn \
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_001/trimmed_A6511_W10_P_R1_001.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_002/trimmed_A6511_W10_P_R1_002.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_003/trimmed_A6511_W10_P_R1_003.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_004/trimmed_A6511_W10_P_R1_004.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_005/trimmed_A6511_W10_P_R1_005.fastq.gz \
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_001/trimmed_A6511_W10_P_R2_001.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_002/trimmed_A6511_W10_P_R2_002.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_003/trimmed_A6511_W10_P_R2_003.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_004/trimmed_A6511_W10_P_R2_004.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_P_005/trimmed_A6511_W10_P_R2_005.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./A6511_W10_P_ \
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
echo "mkdir $HOME/scratch/PPDbRNAseqTimeCourse/STAR2.7.8a_alignment/$foldername; \
cd $HOME/scratch/PPDbRNAseqTimeCourse/STAR2.7.8a_alignment/$foldername; \
STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir /home/ccorreia/scratch/PPDbRNAseqTimeCourse/ARS-UCD1.2/STAR2.7.8a_index \
--readFilesIn $file,$file2,$file3,$file4,$file5 \
$read1,$read2,$read3,$read4,$read5 --readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./${foldername}_ --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done

# Run script on Stampede:
chmod 755 alignment.sh
nohup ./alignment.sh > alignment.sh.nohup &

# Check nohup.out file to see how many jobs finished successfully:
grep -c 'finished successfully' alignment.sh.nohup

# Merge all STAR log.final.out files into a single file:
for file in `find /home/ccorreia/scratch/PPDbRNAseqTimeCourse/STAR2.7.8a_alignment \
-name *Log.final.out`; \
do perl /home/nnalpas/SVN/star_report_opener.pl -report $file; done;

# Transfer All_star_log_final_out.txt file to laptop via SCP.

####################################
# MultiQC for consolidated reports #
####################################

# Required software is MultiQc version 1.9, consult manual/tutorial
# for details: https://multiqc.info/docs/#-20

cd /home/ccorreia/scratch/PPDbRNAseqTimeCourse
multiqc STAR2.7.8a_alignment

#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.9, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment
cd !$

# Create a bash script to perform FastQC quality check on aligned BAM files:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/STAR2.7.8a_alignment \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment $file" >> \
fastqc_aligned.sh; done;

# Split and run all scripts on Stampede
split -d -l 40 fastqc_aligned.sh fastqc_aligned.sh.
for script in `ls fastqc_aligned.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Delete all the HTML files:
rm -r *.html

####################################
# MultiQC for consolidated reports #
####################################

# Required software is MultiQc version 1.9, consult manual/tutorial
# for details: https://multiqc.info/docs/#-20

cd $HOME/scratch/PPDbRNAseqTimeCourse/quality_check
multiqc post_alignment

#########################
# Calculate insert size #
#########################

### Ran Insert size steps on Rodeo, not Stampede.

# Enter working directory:
cd /home/workspace/ccorreia/PPD-b_insert_size

# Sort aligned BAM files:
for file in \
`find /home/workspace/ccorreia/PPD-b_insert_size/STAR2.7.8a_alignment \
-name *Aligned.out.bam`; \
do outfile=`basename $file | perl -p -e 's/Aligned.out.bam/Sorted.out.bam/'`; \
echo "java -jar /usr/local/src/picard/build/libs/picard.jar \
SortSam I=$file O=$outfile SORT_ORDER=coordinate" >> sort.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 40 sort.sh sort.sh.
chmod 755 sort.sh.*
for script in `ls sort.sh.*`;
do nohup ./$script > ${script}.nohup &
done

# Collect insert sizes:
for file in `ls *_Sorted.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Sorted.out.bam//'`; \
echo "java -jar /usr/local/src/picard/build/libs/picard.jar \
CollectInsertSizeMetrics \
I=$file \
O=${sample}_insert_size_metrics.txt \
H=${sample}_insert_size_histogram.pdf M=0.5" >> collect_insert_size.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 40 collect_insert_size.sh collect_insert_size.sh.
chmod 755 collect_insert_size.sh.*
for script in `ls collect_insert_size.sh.*`;
do nohup ./$script > ${script}.nohup &
done

# Collect insert size metrics for all samples into one file:
for file in `ls /home/workspace/ccorreia/PPD-b_insert_size/*_insert_size_metrics.txt`; \
do sample=`basename $file | perl -p -e 's/_insert_size_metrics.txt//'`; \
stats=`sed -n '/MEDIAN_INSERT_SIZE/{n;p;}' $file`; \
printf "${sample}\t${stats}\n" >> All_insert_size.txt; \
done

# Add header to summary stats file:
header=`grep 'MEDIAN_INSERT_SIZE' /home/workspace/ccorreia/PPD-b_insert_size/A6636_W10_U_insert_size_metrics.txt`; \
sed -i $"1 i\Sample_id\t${header}" \
All_insert_size.txt

# Transfer All_insert_size.txt to laptop via SCP.

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################

# Required package is featureCounts, which is part of Subread 1.6.4 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directories:
cd $HOME/scratch/PPDbRNAseqTimeCourse/
mkdir -p Count_summarisation/sense
cd !$

# Run featureCounts with one sample to check if it is working fine:
featureCounts -a \
/home/ccorreia/scratch/PPDbRNAseqTimeCourse/ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.gff \
-B -p -C -R CORE -s 1 -T 15 -t gene -g Dbxref -F GFF -o ./counts.txt \
$HOME/scratch/PPDbRNAseqTimeCourse/STAR2.7.8a_alignment/A6511_W10_P/A6511_W10_P_Aligned.out.bam


# Create a bash script to run featureCounts on BAM file containing multihits and
# uniquely mapped reads using the stranded parameter:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/STAR2.7.8a_alignment \
-name *_Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Aligned.out.bam//'`; \
echo "mkdir $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/$sample; \
cd $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/$sample; \
featureCounts -a \
/home/ccorreia/scratch/PPDbRNAseqTimeCourse/ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.gff \
-B -p -C -R CORE -s 1 -T 15 -t gene -g Dbxref -F GFF \
-o ${sample}_sense-counts.txt $file" >> sense_count.sh; done

# Run script on Stampede:
chmod 755 sense_count.sh
nohup ./sense_count.sh > sense_count.sh.nohup &

# Check if all files were processed:
grep -c 'Summary of counting results can be found in file' sense_count.sh.nohup

# Create bash script to merge stats info from .featureCounts from all samples
# into a single file:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/ \
-name *.featureCounts`; do echo echo \
"\`basename $file\` \`cut $file -f2 | sort | uniq -c | perl -p -e 's/\n/ /'\` >> \
annotation_summary_sense.txt" >> annotation_summary_sense.sh
done

# Run script on Stampede:
chmod 755 annotation_summary_sense.sh
nohup ./annotation_summary_sense.sh &

# Check that all files were processed:
grep -c '.featureCounts' annotation_summary_sense.txt

# Copy all *sense-counts.txt files to temporary folder:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/tmp

for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/ \
-name *sense-counts.txt`; do cp $file \
-t $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/tmp; \
done

# Transfer temporary folder to laptop, then remove it:
rm -r tmp

####################################
# MultiQC for consolidated reports #
####################################

# Required software is MultiQc version 1.9, consult manual/tutorial
# for details: https://multiqc.info/docs/#-20

cd $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation
multiqc sense


