SAMPLES=("oPK1" "oPK2" "oPK3" "oPK4" "oPK5" "yPK1" "yPK2" "yPK3" "yPK4" "yPK5")
SAMPLE="${SAMPLES[${SGE_TASK_ID}-1]}"

#Trimming for bisulfite sequence read file 
fastp -i ${SAMPLE}_Foward.fastq.gz -I ${SAMPLE}_Reverse.fastq.gz -o fastp_${SAMPLE}_Foward.fastq.gz -O fastp_${SAMPLE}_Reverse.fastq.gz -f 10 -F 10

#Merge reference genome and lambda sequence
cat Referance_genome.fa lambda_phage_genome.fa > Reference_lambda.fa

#Map to reference genome
bsmap -r 0 -d Reference_lambda.fa -a fastp_${SAMPLE}_Foward.fastq.gz -b fastp_${SAMPLE}_Reverse.fastq.gz -o ${SAMPLE}_mapped.bam

#Split bam file
bamtools split -in ${SAMPLE}_mapped.bam
bamtools merge -in ${SAMPLE}_mapped.++.bam -in ${SAMPLE}_mapped.+-.bam -out ${SAMPLE}_mapped.top.bam
bamtools merge -in ${SAMPLE}_mapped.-+.bam -in ${SAMPLE}_mapped.--.bam -out ${SAMPLE}_mapped.bottom.bam

#Remove duplicated reads
for strand in top bottom
do
    samtools sort -o ${SAMPLE}_mapped.${strand}.sort.bam ${SAMPLE}_mapped.${strand}.bam
    picard MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=${SAMPLE}_mapped.${strand}.sort.bam \
    OUTPUT=${SAMPLE}_mapped.${strand}.rmdups.bam METRICS_FILE=${SAMPLE}_mapped.${strand}.rmdups.txt \
    REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true
done

bamtools merge -in ${SAMPLE}_mapped.top.rmdups.bam -in ${SAMPLE}_mapped.bottom.rmdups.bam \
-out ${SAMPLE}_mapped.rmdups.bam

#Remove unpaired and unmapped reads
bamtools filter -isMapped true -isPaired true -isProperPair true -forceCompression -in ${SAMPLE}_mapped.rmdups.bam -out ${SAMPLE}_mapped.filter.bam

#Remove overlapping reads
bam clipOverlap --stats --in ${SAMPLE}_mapped.filter.bam --out ${SAMPLE}_mapped.clipped.bam

#Calculate methylation ratio
bsmap mathratio.py -d Reference_lambda.fa -m 1 -o ${SAMPLE}.methylation_results.txt \ ${SAMPLE}_mapped.clipped.bam

bsmap mathratio.py -d Reference_lambda.fa -m 1 -i skip -c lambda \
-o ${SAMPLE}.conversion_ratio.txt ${SAMPLE}_mapped.clipped.bam

for file_name in {o,y}PK{1..5}.conversion_ratio.txt; do; echo $file_name; awk '{a += $7; b += $8} END \
{print (1 - a/b)*100}' $file_name; done

#Calculate total methylation ratio
sed -e "/lambda/d" ${SAMPLE}.methylation_results.txt | awk '{if($4="CG"){{a += $7; b += $8} END \
{print a, b, a/b*100}}}' > ${SAMPLE}.gDNA_CG_total_met_ratio.txt

sed -e "/lambda/d" ${SAMPLE}.methylation_results.txt | awk '{{a += $7; b += $8} END \
{print a, b, a/b*100}}' > ${SAMPLE}.gDNA_total_met_ratio.txt

#Split methylated cytocine
for context in "CG" "CHG" "CHH"
do
    grep "${context}" ${SAMPLE}.mathylation_results.txt > ${SAMPLE}.${context}_meth.txt
done