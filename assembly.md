### Steps for generating short-read genome assembly

1. Read cleaning (fastq)
2. Removing contaminants from reads (kraken)
3. Merging read pairs (bbmerge)  
4. Genome short-read assembly (spades)  


### 1. Read cleaning (if not done before)
fastp v0.24.0

Command for PE reads:
```
fastp --thread ${THREADS} \
-i $PE1 -I $PE2 \
-o 02_trimmed/${SAMPLE}_R1P.fastq.gz -O 02_trimmed/${SAMPLE}_R2P.fastq.gz \
-h 02_trimmed/${SAMPLE}_report.html -j 02_trimmed/${SAMPLE}_report.json \
-n 5 -q 15 -u 40 -l ${MINLEN} \
--detect_adapter_for_pe \
-5 --cut_front_window_size 1 --cut_front_mean_quality 3 \
-r --cut_right_window_size 4 --cut_right_mean_quality 15 \
-g \
--dont_eval_duplication
```

Command for SE reads:
```
fastp --thread ${THREADS} \
-i $SE \
-o 02_trimmed/${SAMPLE}_SE.fastq.gz \
-h 02_trimmed/${SAMPLE}_report.html -j 02_trimmed/${SAMPLE}_report.json \
-n 5 -q 15 -u 40 -l ${MINLEN} \
-5 --cut_front_window_size 1 --cut_front_mean_quality 3 \
-r --cut_right_window_size 4 --cut_right_mean_quality 15 \
-g \
--dont_eval_duplication
```

Specific options:  
MINLEN=31  
Important to use option "g" to trim polyG tails.  

### 2. Removing contaminants
Kraken v2 binaries downloaded locally.  
Downloading Kraken database:
```
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz
```
Command for filtering:
```
${KRAKEN_DIR}/kraken2 --db ${KRAKEN_DB} --confidence 0.1 --threads ${KRAKEN_THREADS} --quick --gzip-compressed --paired \
--unclassified-out ./${DIR}/${SAMPLE}#.fq $PE1 $PE2 \
--output ./${DIR}/${SAMPLE}_out.txt \
--report ./${DIR}/${SAMPLE}_report.txt
```

Specific options:  
confidence=0.1 - if confidence is not set, many reads will be classified as contaminants  

### 3. Merging read pairs
bbmap-39.13-1  
```
bbmerge.sh in1=./02_kraken/${SAMPLE}_1.fq in2=./02_kraken/${SAMPLE}_2.fq out=./${DIR}/${SAMPLE}.fq.gz outu1=./${DIR}/${SAMPLE}_1.fq.gz outu2=./${DIR}/${SAMPLE}_2.fq.gz
```

### 4. Genome assembly
Spades 4.0.0  
```
spades.py --isolate -1 ./03_bbmerge/${SAMPLE}_1.fq.gz -2 ./03_bbmerge/${SAMPLE}_2.fq.gz --merged ./03_bbmerge/${SAMPLE}.fq.gz -o ./${DIR}/spades_${SAMPLE} -t ${SPADES_THREADS} -m ${SPADES_MEMORY_LIMIT}
```
