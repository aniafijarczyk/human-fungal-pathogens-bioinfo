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

Specific options used:  
MINLEN=31  
Important to use option "g" to trim polyG tails.  

### 2. Removing contaminants

Kraken v2 binaries downloaded locally.  
Downloading Kraken database:
```
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz
```
Command for filtering
