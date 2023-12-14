"""
For trimming with TRIMMOMATIC, trimmomatic.jar must be in directory folder

STAR was used for alingment and the reference index was generated using :
STAR --runMode genomeGenerate \
--runThreadN 6 \
--genomeDir ./ \
--genomeFastaFiles ./CHIKFLIC_genome.fa" \
--genomeSAindexNbases 6


"""
# loop in file folder to find R1 and R2 fastq files and run Trimmomatic
for R1 in *_R1_001.fastq.gz; do
R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

OUTPUT_R1_PAIRED="trimmed_${R1/_R1_001.fastq.gz/_R1_001_paired.fastq.gz}"
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
