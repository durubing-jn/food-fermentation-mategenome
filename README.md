# food-fermentation-mategenome

## 1. Raw data download

```sh
fastq-dump ACCESSION ID --split-files
```

## 2. Quality control

```SH
#pair-end metageomic sequencing
trim_galore -o OUR_FOLDER --gzip --paired RAW.DATA_1.fastq  RAW.DATA_2.fastq
#single-end metageomic sequencing
trim_galore -o OUR_FOLDER --gzip RAW.DATA.fastq
#quality check
fastqc -t 24 -o OUR_FOLDER -f fastq CLEAN.DATA.fastq
```

