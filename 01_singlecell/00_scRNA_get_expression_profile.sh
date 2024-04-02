###get expression profile from scRNA-seq fastq file

#!/bin/bash

#slurm options
#SBATCH -p intel-debug                   
#SBATCH -q debug                        
#SBATCH -J model                         
#SBATCH -c 10                                  
#SBATCH -o model.log                    

id=$id
sample_name=$sample_name
fqdir=$fqdir
path=$path
result=$path/result

read1=$id.read1
read2=$id.read2
read2_newbarcode_fq=$read2.fq
clean_data_read1=$fqdir/$sample_name/${sample_name}_1.clean.fq.gz
clean_data_read2=$fqdir/$sample_name/${sample_name}_2.clean.fq.gz

source ~/.bashrc
module load umi_tools/1.1.1
module load star/2.7.8a
module load samtools/1.13
module load subread/2.0.2
module load fastqc/0.11.9

umi_tools=$path_to_umi_tools
STAR=$path_to_STAR
samtools=$path_to_samtools
featureCounts=$path_to_featureCounts
trim_galore=$path_to_trim_galore
fastqc=$path_to_fastqc
# reference
mm10gtf=$path_to_gtf
mm10index=$path_to_starindex

if [ ! -d "$result" ]; then
  mkdir $result
fi
if [ ! -d "$result/01_filter" ]; then
  mkdir $result/01_filter
fi
if [ ! -d "$result/01_fastqc" ]; then
  mkdir $result/01_fastqc
fi
if [ ! -d "$result/01_fastqc2" ];then
 mkdir $result/01_fastqc2
fi
if [ ! -d "$result/02_add_barcode" ]; then
  mkdir $result/02_add_barcode
fi
if [ ! -d "$result/03_extract_barcode" ]; then
  mkdir $result/03_extract_barcode
fi

if [ ! -d "$result/04_mapping_STAR" ]; then
  mkdir $result/04_mapping_STAR
fi

if [ ! -d "$result/05_featureCount" ]; then
  mkdir $result/05_featureCount
fi

if [ ! -d "$result/06_umiCount" ]; then
  mkdir $result/06_umiCount
fi

#fastqc
$fastqc $clean_data_read1 -d $result/01_fastqc -t 10 --noextract -o $result/01_fastqc
$fastqc $clean_data_read2 -d $result/01_fastqc -t 10 --noextract -o $result/01_fastqc

# # filter q20
$trim_galore --paired -q 20 --fastqc -o $result/01_filter $clean_data_read1 $clean_data_read2

#fastqc
id1=$result/01_filter/${sample_name}_1.clean_val_1.fq.gz
id2=$result/01_filter/${sample_name}_2.clean_val_2.fq.gz

$fastqc $id1 -d $result/01_fastqc2 -t 10 --noextract -o $result/01_fastqc2
$fastqc $id2 -d $result/01_fastqc2 -t 10 --noextract -o $result/01_fastqc2

gunzip -c $result/01_filter/${sample_name}_1.clean_val_1.fq.gz > $result/01_filter/$read1.fq
gunzip -c $result/01_filter/${sample_name}_2.clean_val_2.fq.gz > $result/01_filter/$read2.fq

# # make new barcode
sed '2~2s/^/AA/g' $result/01_filter/$read2.fq > $result/02_add_barcode/$read2_newbarcode_fq
# extract barcode
stdin_file=$result/02_add_barcode/$read2_newbarcode_fq
stdot_file=$result/03_extract_barcode/$read2.extract_barcodes.fq
read2_in_file=$result/01_filter/$read1.fq
read2_out_file=$result/03_extract_barcode/$read1.extract_barcodes.fq
whitelist_file=$result/03_extract_barcode/whitelist_$id.txt

umi_tools whitelist  --stdin $stdin_file  --bc-pattern=CCCCCCCCNNNNNNNN  --set-cell-number=110  --log2stderr > $whitelist_file
umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN  --stdin $stdin_file  --stdout $stdot_file  --read2-in $read2_in_file  --read2-out=$read2_out_file  --filter-cell-barcode --whitelist=$whitelist_file

# mapping
$STAR --runThreadN 10 --genomeDir $mm10index --readFilesIn $read2_out_file --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix $result/04_mapping_STAR/$id.

# featurecount
featurecount_outfile=$result/05_featureCount/$id.gene_assigned
featurecount_outbam=$result/04_mapping_STAR/$id.Aligned.sortedByCoord.out.bam
featureCounts -a $mm10gtf -o $featurecount_outfile -R BAM $featurecount_outbam      

# resort bam file
input_file=$result/05_featureCount/$id.Aligned.sortedByCoord.out.bam.featureCounts.bam
sorted_file=$result/05_featureCount/$id.assigned_sorted.bam
samtools sort $input_file -o $sorted_file
samtools index $sorted_file

# umicount
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I $sorted_file -S $result/06_umiCount/$id.counts.tsv.gz --wide-format-cell-counts

