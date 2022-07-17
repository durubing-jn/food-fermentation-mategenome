# Food fermentation mategenome

## 1. Raw data download

```sh
for i in $(cat ACCESSION_ID); do
fastq-dump ${i} --split-files
done
```

## 2. Quality control

```SH
#pair-end metageomic sequencing
trim_galore -o OUTPUT_QC_FOLDER --gzip --paired RAW.DATA_1.fastq  RAW.DATA_2.fastq

#single-end metageomic sequencing
trim_galore -o OUTPUT_QC_FOLDER --gzip RAW.DATA.fastq

#quality check
fastqc -t 24 -o OUR_FOLDER -f fastq CLEAN.DATA_1.fastq.gz
fastqc -t 24 -o OUR_FOLDER -f fastq CLEAN.DATA.fastq.gz
```

### 3. Assembly

```sh
#MEGAHIT, pair-end metageomic sequencing
megahit -1 CLEAN.DATA_1.fastq.gz -2 CLEAN.DATA_1.fastq.gz -t 30  -o OUTPUT_ASSEMBLY_MEGAHIT_PAIR

#MEGAHIT, single-end metageomic sequencing
megahit -r CLEAN.DATA.fastq.gz -t 30  -o OUTPUT_ASSEMBLY_MEGAHIT_SINGLE

#metaspades, pair-end metageomic sequencing
metaspades.py --tmp-dir TEM_OUTPUT_FOLDER -t 30 -o OUTPUT_ASSEMBLY_metaspades_PAIR -1 CLEAN.DATA_1.fastq.gz -2 CLEAN.DATA_1.fastq.gz

#metaspades, single-end metageomic sequencing
metaspades.py --tmp-dir 2.2_assembly-metaspades/tmp/${i}.tmp -t 4 -o OUTPUT_ASSEMBLY_metaspades_SINGLE --s1 CLEAN.DATA.fastq.gz
```

## 4. Binning

```sh
#Filter contigs by length (1500 bp)
seqtk seq -L 1500 3.2_assembly_results-metaspades/${i}_scaffolds.fasta > 4.2_assembly_file_results_metaspades/${i}.f_fa

#Binning analysis by metabat2
for i in DRR205701 DRR205702 DRR205703; do
mkdir ${i}
bowtie2-build -f ~/1_fermentated_food/4.3_co-assembly_file_results/T_Chineseliquor.fa.f ${i}/${i}_final --threads 10
bowtie2 -1 ~/1_fermentated_food/1_clean_data/${i}_1_val_1.fq.gz -2 ~/1_fermentated_food/1_clean_data/${i}_2_val_2.fq.gz -p 10 -x ${i}/${i}_final -S ${i}/${i}_final.sam 2>0_SAM/${i}.bowtie2.stat
samtools view -@ 10 -b -S ${i}/${i}_final.sam -o ${i}/${i}_final.bam
samtools sort -@ 10 -l 9 -O BAM ${i}/${i}_final.bam -o ${i}/${i}_final.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth ${i}/${i}_final.depth.txt ${i}/${i}_final.sorted.bam
metabat2 -m 1500 -t 10 -i ~/1_fermentated_food/4.3_co-assembly_file_results/T_Chineseliquor.fa.f -a ${i}/${i}_final.depth.txt -o 1_bins/${i} -v
done

#Bing analysis by metabin2
for i in DRR205701 DRR205702 DRR205703; do
genomeCoverageBed -ibam ${i}/${i}_final.sorted.bam > ${i}/$i.histogram.tab
python ~/1_fermentated_food/calculate-contig-coverage.py ${i}/$i.histogram.tab
run_MaxBin.pl -contig ~/1_fermentated_food/4.3_co-assembly_file_results/T_Chineseliquor.fa.f -abund ${i}/${i}.histogram.tab.coverage.tab -max_iteration 50 -out 1.1_bin-metabin2/${i}_metabin2 -thread 10
done

#Bing analysis by concoct
for i in DRR205702 DRR205703; do
cut_up_fasta.py ~/1_fermentated_food/4.3_co-assembly_file_results/T_Chineseliquor.fa.f -c 10000 -o 0 --merge_last -b ${i}/${i}.contigs_10K.bed > ${i}/${i}.contigs_10K.fa
samtools index 2.2_final.sorted.bam_metahit_co_assembly/${i}_final.sorted.bam -@ 12
concoct_coverage_table.py ${i}/${i}.contigs_10K.bed 2.2_final.sorted.bam_metahit_co_assembly/${i}_final.sorted.bam > ${i}/${i}.coverage_table.tsv
concoct --composition_file ${i}/${i}.contigs_10K.fa --coverage_file ${i}/${i}.coverage_table.tsv -b ${i}/${i}.concoct_output --threads 12
merge_cutup_clustering.py ${i}/${i}.concoct_output_clustering_gt1000.csv > ${i}/${i}.clustering_merged.csv
extract_fasta_bins.py ~/1_fermentated_food/4.3_co-assembly_file_results/T_Chineseliquor.fa.f ${i}/${i}.clustering_merged.csv --output_path 1.3_bin-concoct/${i}
done
```

## 5. DAS_TOOL

```SH
Fasta_to_Scaffolds2Bin.sh -i 1.1.2_bins-metabat2/${i} -e fa > 2.das_tool/${i}_metabat2.scaffolds2bin.tsv
DAS_Tool -i 2.das_tool/${i}_metabat2.scaffolds2bin.tsv,2.das_tool/${i}_metabin2.scaffolds2bin.tsv,2.das_tool/${i}_concoct.scaffolds2bin.tsv -l metabat,maxbin,concoct -c ~/1_fermentated_food/4_assembly_file_results/${i}.f_fa -o 3.3das_tool_results/${i} --search_engine diamond --write_bins 1 --score_threshold 0 -t 30
```

## 6.CheckM

```sh
checkm lineage_wf -x fa 0.genome 1.checkm -t 30 --tmpdir bin_checkm.tmp
```

## 7.dRep

```sh
dRep dereplicate 2.dRep/drep_${i}/ -g 15_fermented_foods/${i}/*.fa -p 2 --ignoreGenomeQuality -pa 0.95 -sa 0.99 --S_algorithm fastANI

dRep dereplicate -g 2.dRep/drep_${i}/dereplicated_genomes/*.fa -p 20 3.SGB/${i}/ --ignoreGenomeQuality -pa 0.9 -sa 0.95 --S_algorithm fastANI
```

## 8. GTDB_tk

```sh
gtdbtk classify_wf --genome_dir 01.nr_99 --out_dir 05.genome_gtdbtk --cpus 30 --scratch_dir gtdbtk.temp --pplacer_cpus 32 --extension fa
```

## 9. Phylogenetic analysis

```sh
bmge -i gtdbtk.bac120.user_msa.fasta -t AA -g 0.5 -h 1 -b 1 -w 1 -of out_bmge_trimmed_2.fasta
FastTree out_bmge_trimmed_2.fasta > Bacteria_trimmed_2.tree
```

## 10.  antismash analysis

```sh
antismash 03.nr_genome/01.nr_99/MAG_rename/ACT/${1}.fa --taxon bacteria --output-dir 05.SM_BGC/01.antismash_fullhmmer/${1} --genefinding-tool prodigal --cb-knownclusters -c 3 --cc-mibig --fullhmmer
```

## 11. Bio-SCAPE

```sh
python3 ~/1_fermentated_food/7.genome_analysis/5.SM_ANALYSIS/BiG-SCAPE-master/bigscape.py -i input -o output_20220615 --pfam_dir ~/1_fermentated_food/7.genome_analysis/5.SM_ANALYSIS/BiG-SCAPE-master/hmm -c 30 --cutoffs 0.7 --include_singletons --mode auto --mix --no_classify --mibig --verbose > run.log
```

## 12. Bio-SLiCE

```sh
bigslice -i 02.cleandata 03.humangut.bigslice
bigslice --query ~/1_fermentated_food/7.genome_analysis/05.SM_BGC/02.BiG-SCAPE/input --n_ranks 1 full_run_result -t 30
```

## 13. BGC_function_prediction

```sh
./cluster_function_prediction.py BGCs_for_prediction/antiSMASH/MAG_1_228.region001.gbk BGCs_for_prediction/CARD/MAG_1.txt --output test --antismash_version 5 --rgi_version 5
```

