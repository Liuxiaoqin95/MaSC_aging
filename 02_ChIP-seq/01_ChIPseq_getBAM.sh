###get mapped bam file from ChIPseq fastq file

#!/bin/bash

#slurm options
#SBATCH -p intel-debug                   
#SBATCH -q debug                               
#SBATCH -J chipseq                         
#SBATCH -c 10                                   
#SBATCH -o chipseq.log                        


dir_fq=$fastq_dir
id0=$result_file_id
id1=$read1
id2=$read2
index=$bowtie2_index  
result=$result_dir

bwa=$path_to_bwa
samtools=$path_to_samtools
java=$path_to_java
MarkDuplicates=$path_to_MarkDuplicates
fastqc=$path_to_fastqc
trim_galore=$path_to_trim_galore
bowtie2=$path_to_bowtie2
fq1=$dir_fq/$id1
fq2=$dir_fq/$id2
pre=${id1/_R1_001.fastq.gz/}
post1="_R1_001_val_1.fq.gz"
post2="_R2_001_val_2.fq.gz"

$trim_galore --path_to_cutadapt=~/.local/bin/cutadapt --paired -q 20 --fastqc -o $result/01_Filter_q20 $fq1 $fq2 

$fastqc $result/01_Filter_q20/${pre}${post1} -d $result/01_fastqc -t 12 --noextract -o $result/01_fastqc
$fastqc $result/01_Filter_q20/${pre}${post2} -d $result/01_fastqc -t 12 --noextract -o $result/01_fastqc

$bowtie2 -x $index --no-1mm-upfront -1 $result/01_Filter_q20/${pre}${post1} -2 $result/01_Filter_q20/${pre}${post2} -p 30 -S $result/02_bowtie2Mapping/${id0}.raw.sam

## find bad cigar read names 
cat $result/02_bowtie2Mapping/${id0}.raw.sam | awk 'BEGIN {FS="\t"; OFS="\t"} !/^@/ && $6!="*" {cigar=$6; gsub("[0-9]+D","",cigar); n=split(cigar,vals,"[A-Z]"); s=0; for(i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10); if(s!=seqlen)print $1"\t";}' | sort | uniq > ${id0}.badReads.tmp

if [[ $(cat ${id0}.badReads.tmp|wc -l) -gt 0 ]]
then
	cat $result/02_bowtie2Mapping/${id0}.raw.sam | grep -v -F -f ${id0}.badReads.tmp | $samtools view -Su - | $samtools sort -o $result/02_bowtie2Mapping/${id0}.raw.bam
else
	$samtools view -Su $result/02_bowtie2Mapping/${id0}.raw.sam | $samtools sort -o $result/02_bowtie2Mapping/${id0}.raw.bam
fi

$samtools index $result/02_bowtie2Mapping/${id0}.raw.bam $result/02_bowtie2Mapping/${id0}.raw.bai
$samtools flagstat $result/02_bowtie2Mapping/${id0}.raw.bam > $result/02_bowtie2Mapping/${id0}.raw.flagstat.qc
## 01. alignments end

## 02. post alignments start
$samtools view -F 1804 -f 2 -q 30 -u $result/02_bowtie2Mapping/${id0}.raw.bam |$samtools sort -n -o $result/02_bowtie2Mapping/tmp.${id0}.filt.srt.nmsrt
$samtools fixmate -r $result/02_bowtie2Mapping/tmp.${id0}.filt.srt.nmsrt - | $samtools view -F 1804 -f 2 -h - -O BAM| $samtools sort -o $result/02_bowtie2Mapping/${id0}.filt.srt.bam
rm $result/02_bowtie2Mapping/tmp.${id0}.filt.srt.nmsrt

$java -jar $MarkDuplicates INPUT=$result/02_bowtie2Mapping/${id0}.filt.srt.bam OUTPUT=$result/02_bowtie2Mapping/${id0}.filt.srt.rmdup.bam METRICS_FILE=$result/02_bowtie2Mapping/${id0}.rmdup.metrics VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
mv $result/02_bowtie2Mapping/${id0}.filt.srt.rmdup.bam $result/02_bowtie2Mapping/${id0}.filt.srt.bam
$samtools index $result/02_bowtie2Mapping/${id0}.filt.srt.bam $result/02_bowtie2Mapping/${id0}.filt.srt.bai
$samtools view -F 1804 -f 2 -b $result/02_bowtie2Mapping/${id0}.filt.srt.bam > $result/03_final/${id0}.filt.srt.nodup.bam
$samtools sort -n -o $result/03_final/${id0}.filt.nmsrt.nodup.bam $result/03_final/${id0}.filt.srt.nodup.bam 
$samtools index $result/03_final/${id0}.filt.srt.nodup.bam $result/03_final/${id0}.filt.srt.nodup.bai
$samtools flagstat $result/03_final/${id0}.filt.srt.nodup.bam > $result/03_final/${id0}.filt.srt.nodup.flagstat.qc


rm $result/02_bowtie2Mapping/${id0}.raw.bam $result/02_bowtie2Mapping/${id0}.raw.bai $result/02_bowtie2Mapping/${id0}.rmdup.metrics $result/02_bowtie2Mapping/${id0}.raw.flagstat.qc
rm $result/02_bowtie2Mapping/${id0}.raw.sam
rm ${id0}.badReads.tmp

