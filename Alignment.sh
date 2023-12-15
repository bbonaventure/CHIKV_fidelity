"""
For trimming with TRIMMOMATIC, trimmomatic.jar must be in directory folder

STAR was used for alingment and the reference index was generated using :
STAR --runMode genomeGenerate \
--runThreadN 6 \
--genomeDir ./ \
--genomeFastaFiles ./CHIKFLIC_genome.fa" \
--genomeSAindexNbases 6
"""
For deduplication with picard.jar, picard.jar must be in directory folder
"""
# loop in file folder to find R1 and R2 fastq files and run Trimmomatic
for R1 in *_R1_001.fastq.gz; do
R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

OUTPUT_R1_PAIRED="trimmed_${R1/_R1_001.fastq.gz/_R1_001_paired.fastq .gz}"
OUTPUT_R1_UNPAIRED="trimmed_${R1/_R1_001.fastq.gz/_R1_001_unpaired.fastq.gz}"
OUTPUT_R2_PAIRED="trimmed_${R2/_R2_001.fastq.gz/_R2_001_paired.fastq.gz}"
OUTPUT_R2_UNPAIRED="trimmed_${R2/_R2_001.fastq.gz/_R2_001_unpaired.fastq.gz}"
#run trimmomatic
java -jar trimmomatic-0.39.jar PE -phred33 \
"$R1" "$R2" \
"$OUTPUT_R1_PAIRED" "$OUTPUT_R1_UNPAIRED" \
"$OUTPUT_R2_PAIRED" "$OUTPUT_R2_UNPAIRED" \
ILLUMINACLIP:"$ADAPTERS":2:30:10:8:TRUE \
LEADING:$LEADING \
TRAILING:$TRAILING \
SLIDINGWINDOW:$SLIDINGWINDOW \
MINLEN:$MINLEN
done


#take trimmed and paired fastq files to run STAR

for R1 in *_R1_001_paired.fastq.gz
do
R2=${R1/_R1_001_paired.fastq.gz/_R2_001_paired.fastq.gz}

OUT_PREFIX="${R1%_R1_001_paired.fastq.gz}"

# Run STAR alignment
STAR --genomeDir ./ \
--readFilesIn "$R1" "$R2" \
--readFilesCommand gunzip -c \
--runThreadN 8 \
--outReadsUnmapped Fastx \
--outFileNamePrefix "$OUT_PREFIX"
done

##Convertion of SAM files to BAM files and sorting.
mkdir -p BAM_files   
for samfile in *.sam; do                 
base_name=$(basename "$samfile" .sam)                                       
samtools view -bS "$samfile" | samtools sort -o "BAM_files/${base_name}_sorted.bam"
done

#Deduplication loop with picard.jar
for input_file in *.bam; do
output_file="deduplicated_${input_file}"
metrics_file="metrics_${input_file}.txt"

java -jar picard.jar MarkDuplicates \
INPUT="${input_file}" \
OUTPUT="${output_file}" \
METRICS_FILE="${metrics_file}" \
REMOVE_DUPLICATES=true
done

##Indexing of sorted, deduplicated BAM files
for file in deduplicated*.bam; do
samtools index "$file"
done

