# Food fermentation mategenome

## 1. Raw data download

```sh
#
for i in $(cat data/ACCESSION_ID ); do
fastq-dump ${i} --split-files
done
```

## 2. Quality control

```SH
#pair-end metageomic sequencing
for i in $(cat data/pair_end_sample_id ); do
trim_galore -o 01.cleandata --gzip --paired 00.rawdata/${i}_1.fastq  00.rawdata/${i}_2.fastq
done

#single-end metageomic sequencing
for i in $(cat data/single_end_sample_id ); do
trim_galore -o 01.cleandata --gzip 00.rawdata/${i}_1.fastq
done

#quality check
for i in $(cat data/pair_end_sample_id ); do
fastqc -t 24 -o 02.QC -f fastq 01.cleandata/${i}_1_val_1.fq.gz
fastqc -t 24 -o 02.QC -f fastq 01.cleandata/${i}_2_val_2.fq.gz
done

for i in $(cat data/single_end_sample_id ); do
fastqc -t 24 -o 02.QC -f fastq 01.cleandata/${i}_trimmed_1.fastq.gz
done
```

### 3. Assembly

```sh
#MEGAHIT, pair-end metageomic sequencing
for in $(cat data/pair_end_sample_id ); do
megahit -1 01.cleandata/${i}_1_val_1.fq.gz -2 01.cleandata/${i}_2_val_2.fq.gz -t 30  -o 03.assmebly/00.assmebly_megahit/${i}
done

#MEGAHIT, single-end metageomic sequencing
for i in $(cat data/single_end_sample_id ); do
megahit -r 01.cleandata/${i}_trimmed_1.fastq.gz -t 30  -o 03.assmebly/00.assmebly_megahit/${i}

#metaspades, pair-end metageomic sequencing
for in $(cat data/pair_end_sample_id ); do
metaspades.py --tmp-dir TEM_OUTPUT_FOLDER/${i}.tmp -t 30 -o 03.assmebly/01.assmebly_metspades/${i} -1 01.cleandata/${i}_1_val_1.fq.gz -2 01.cleandata/${i}_2_val_2.fq.gz
done

#metaspades, single-end metageomic sequencing
for i in $(cat data/single_end_sample_id); do
metaspades.py --tmp-dir TEM_OUTPUT_FOLDER//${i}.tmp -t 30 -o 03.assmebly/01.assmebly_metspades/${i} --s1 01.cleandata/${i}_trimmed_1.fastq.gz
done
```

## 4. Binning

```sh
#Filter contigs by length (1500 bp)
for i in $(cat data/all_sampe_id); do
seqtk seq -L 1500 03.assmebly/00.assmebly_megahit/${i}/final.contigs.fa > 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa
done

for i in $(cat data/all_sampe_id); do
seqtk seq -L 1500 03.assmebly/01.assmebly_metspades/${i}/contigs.fasta > 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa
done

####################Binning analysis based on metahit assmebly#############################
#pair-end metageomic sequencing
for in $(cat data/pair_end_sample_id ); do
mkdir 05.binning/metahit_bin/${i}
bowtie2-build -f 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa --threads 10
bowtie2 -1 01.cleandata/${i}_1_val_1.fq.gz -2 01.cleandata/${i}_2_val_2.fq.gz -p 10 -x 05.binning/metahit_bin/${i}/${i}_final -S 05.binning/metahit_bin/${i}/${i}_final.sam 2>0_SAM/${i}.bowtie2.stat
samtools view -@ 10 -b -S 05.binning/metahit_bin/${i}/${i}_final.sam -o 05.binning/metahit_bin/${i}/${i}_final.bam
samtools sort -@ 10 -l 9 -O BAM 05.binning/metahit_bin/${i}/${i}_final.bam -o 05.binning/metahit_bin/${i}/${i}_final.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth 05.binning/metahit_bin/${i}/${i}_final.depth.txt 05.binning/metahit_bin/${i}/${i}_final.sorted.bam
metabat2 -m 1500 -t 10 -i 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa -a 05.binning/metahit_bin/${i}/${i}_final.depth.txt -o 05.binning/results/00.metabat2.metahit.bins/${i} -v

genomeCoverageBed -ibam 05.binning/metahit_bin/${i}/${i}_final.sorted.bam > 05.binning/metahit_bin/${i}/${i}.histogram.tab
python ./calculate-contig-coverage.py 05.binning/metahit_bin/${i}/$i.histogram.tab
run_MaxBin.pl -contig 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa -abund 05.binning/metahit_bin/${i}/${i}.histogram.tab.coverage.tab -max_iteration 50 -out 05.binning/results/01.metabin2.metahit.bins/${i}_metabin2 -thread 10

cut_up_fasta.py 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa -c 10000 -o 0 --merge_last -b 05.binning/metahit_bin/${i}/${i}.contigs_10K.bed > 05.binning/metahit_bin/${i}/${i}.contigs_10K.fa
samtools index 05.binning/metahit_bin/${i}/${i}_final.sorted.bam -@ 10
concoct_coverage_table.py 05.binning/metahit_bin/${i}/${i}.contigs_10K.bed 05.binning/metahit_bin/${i}/${i}_final.sorted.bam > 05.binning/metahit_bin/${i}/${i}.coverage_table.tsv
concoct --composition_file 05.binning/metahit_bin/${i}/${i}.contigs_10K.fa --coverage_file 05.binning/metahit_bin/${i}/${i}.coverage_table.tsv -b 05.binning/metahit_bin/${i}/${i}.concoct_output --threads 10
merge_cutup_clustering.py 05.binning/metahit_bin/${i}/${i}.concoct_output_clustering_gt1000.csv > 05.binning/metahit_bin/${i}/${i}.clustering_merged.csv
extract_fasta_bins.py 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa 05.binning/metahit_bin/${i}/${i}.clustering_merged.csv --output_path 05.binning/results/02.concoct.metahit.bins/${i}
done

#single-end metageomic sequencing
for in $(cat data/single_end_sample_id ); do
mkdir 05.binning/metahit_bin/${i}
bowtie2-build -f 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa --threads 10
bowtie2 -U 01.cleandata/${i}_trimmed_1.fastq.gz -p 10 -x 05.binning/metahit_bin/${i}/${i}_final -S 05.binning/metahit_bin/${i}/${i}_final.sam 2>0_SAM/${i}.bowtie2.stat
samtools view -@ 10 -b -S 05.binning/metahit_bin/${i}/${i}_final.sam -o 05.binning/metahit_bin/${i}/${i}_final.bam
samtools sort -@ 10 -l 9 -O BAM 05.binning/metahit_bin/${i}/${i}_final.bam -o 05.binning/metahit_bin/${i}/${i}_final.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth 05.binning/metahit_bin/${i}/${i}_final.depth.txt 05.binning/metahit_bin/${i}/${i}_final.sorted.bam
metabat2 -m 1500 -t 10 -i 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa -a 05.binning/metahit_bin/${i}/${i}_final.depth.txt -o 05.binning/results/00.metabat2.metahit.bins/${i} -v

genomeCoverageBed -ibam 05.binning/metahit_bin/${i}/${i}_final.sorted.bam > 05.binning/metahit_bin/${i}/${i}.histogram.tab
python ./calculate-contig-coverage.py 05.binning/metahit_bin/${i}/$i.histogram.tab
run_MaxBin.pl -contig 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa -abund 05.binning/metahit_bin/${i}/${i}.histogram.tab.coverage.tab -max_iteration 50 -out 05.binning/results/01.metabin2.metahit.bins/${i}_metabin2 -thread 10

cut_up_fasta.py 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa -c 10000 -o 0 --merge_last -b 05.binning/metahit_bin/${i}/${i}.contigs_10K.bed > 05.binning/metahit_bin/${i}/${i}.contigs_10K.fa
samtools index 05.binning/metahit_bin/${i}/${i}_final.sorted.bam -@ 10
concoct_coverage_table.py 05.binning/metahit_bin/${i}/${i}.contigs_10K.bed 05.binning/metahit_bin/${i}/${i}_final.sorted.bam > 05.binning/metahit_bin/${i}/${i}.coverage_table.tsv
concoct --composition_file 05.binning/metahit_bin/${i}/${i}.contigs_10K.fa --coverage_file 05.binning/metahit_bin/${i}/${i}.coverage_table.tsv -b 05.binning/metahit_bin/${i}/${i}.concoct_output --threads 10
merge_cutup_clustering.py 05.binning/metahit_bin/${i}/${i}.concoct_output_clustering_gt1000.csv > 05.binning/metahit_bin/${i}/${i}.clustering_merged.csv
extract_fasta_bins.py 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa 05.binning/metahit_bin/${i}/${i}.clustering_merged.csv --output_path 05.binning/results/02.concoct.metahit.bins/${i}
done


####################Binning analysis based on metaspades assmebly#############################
#pair-end metageomic sequencing
for in $(cat data/pair_end_sample_id ); do
mkdir 05.binning/metaspades_bin/${i}
bowtie2-build -f 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa --threads 10
bowtie2 -1 01.cleandata/${i}_1_val_1.fq.gz -2 01.cleandata/${i}_2_val_2.fq.gz -p 10 -x 05.binning/metaspades_bin/${i}/${i}_final -S 05.binning/metaspades_bin/${i}/${i}_final.sam 2>0_SAM/${i}.bowtie2.stat
samtools view -@ 10 -b -S 05.binning/metaspades_bin/${i}/${i}_final.sam -o 05.binning/metaspades_bin/${i}/${i}_final.bam
samtools sort -@ 10 -l 9 -O BAM 05.binning/metaspades_bin/${i}/${i}_final.bam -o 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth 05.binning/metaspades_bin/${i}/${i}_final.depth.txt 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam
metabat2 -m 1500 -t 10 -i 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa -a 05.binning/metaspades_bin/${i}/${i}_final.depth.txt -o 05.binning/results/00.metabat2.metaspades.bins/${i} -v

genomeCoverageBed -ibam 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam > 05.binning/metaspades_bin/${i}/${i}.histogram.tab
python ./calculate-contig-coverage.py 05.binning/metaspades_bin/${i}/$i.histogram.tab
run_MaxBin.pl -contig 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa -abund 05.binning/metaspades_bin/${i}/${i}.histogram.tab.coverage.tab -max_iteration 50 -out 05.binning/results/01.metabin2.metaspades.bins/${i}_metabin2 -thread 10

cut_up_fasta.py 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa -c 10000 -o 0 --merge_last -b 05.binning/metaspades_bin/${i}/${i}.contigs_10K.bed > 05.binning/metaspades_bin/${i}/${i}.contigs_10K.fa
samtools index 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam -@ 10
concoct_coverage_table.py 05.binning/metaspades_bin/${i}/${i}.contigs_10K.bed 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam > 05.binning/metaspades_bin/${i}/${i}.coverage_table.tsv
concoct --composition_file 05.binning/metaspades_bin/${i}/${i}.contigs_10K.fa --coverage_file 05.binning/metaspades_bin/${i}/${i}.coverage_table.tsv -b 05.binning/metaspades_bin/${i}/${i}.concoct_output --threads 10
merge_cutup_clustering.py 05.binning/metaspades_bin/${i}/${i}.concoct_output_clustering_gt1000.csv > 05.binning/metaspades_bin/${i}/${i}.clustering_merged.csv
extract_fasta_bins.py 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa 05.binning/metaspades_bin/${i}/${i}.clustering_merged.csv --output_path 05.binning/results/02.concoct.metaspades.bins/${i}
done

#single-end metageomic sequencing
for in $(cat data/single_end_sample_id ); do
mkdir 05.binning/metaspades_bin/${i}
bowtie2-build -f 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa --threads 10
bowtie2 -U 01.cleandata/${i}_trimmed_1.fastq.gz -p 10 -x 05.binning/metaspades_bin/${i}/${i}_final -S 05.binning/metaspades_bin/${i}/${i}_final.sam 2>0_SAM/${i}.bowtie2.stat
samtools view -@ 10 -b -S 05.binning/metaspades_bin/${i}/${i}_final.sam -o 05.binning/metaspades_bin/${i}/${i}_final.bam
samtools sort -@ 10 -l 9 -O BAM 05.binning/metaspades_bin/${i}/${i}_final.bam -o 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth 05.binning/metaspades_bin/${i}/${i}_final.depth.txt 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam
metabat2 -m 1500 -t 10 -i 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa -a 05.binning/metaspades_bin/${i}/${i}_final.depth.txt -o 05.binning/results/00.metabat2.metaspades.bins/${i} -v

genomeCoverageBed -ibam 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam > 05.binning/metaspades_bin/${i}/${i}.histogram.tab
python ./calculate-contig-coverage.py 05.binning/metaspades_bin/${i}/$i.histogram.tab
run_MaxBin.pl -contig 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa -abund 05.binning/metaspades_bin/${i}/${i}.histogram.tab.coverage.tab -max_iteration 50 -out 05.binning/results/01.metabin2.metaspades.bins/${i}_metabin2 -thread 10

cut_up_fasta.py 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa -c 10000 -o 0 --merge_last -b 05.binning/metaspades_bin/${i}/${i}.contigs_10K.bed > 05.binning/metaspades_bin/${i}/${i}.contigs_10K.fa
samtools index 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam -@ 10
concoct_coverage_table.py 05.binning/metaspades_bin/${i}/${i}.contigs_10K.bed 05.binning/metaspades_bin/${i}/${i}_final.sorted.bam > 05.binning/metaspades_bin/${i}/${i}.coverage_table.tsv
concoct --composition_file 05.binning/metaspades_bin/${i}/${i}.contigs_10K.fa --coverage_file 05.binning/metaspades_bin/${i}/${i}.coverage_table.tsv -b 05.binning/metaspades_bin/${i}/${i}.concoct_output --threads 10
merge_cutup_clustering.py 05.binning/metaspades_bin/${i}/${i}.concoct_output_clustering_gt1000.csv > 05.binning/metaspades_bin/${i}/${i}.clustering_merged.csv
extract_fasta_bins.py 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa 05.binning/metaspades_bin/${i}/${i}.clustering_merged.csv --output_path 05.binning/results/02.concoct.metaspades.bins/${i}
done
```

## 5. DAS_TOOL

```SH
#Preparing configuration Files
for i in $(cat data/all_sampe_id); do
cp 05.binning/results/00.metabat2.metahit.bins/${i}*.fa 05.binning/bins_metahit/${i}_metabat2/
cp 05.binning/results/01.metabin2.metahit.bins/${i}_metabin2.*.fasta 05.binning/bins_metahit/${i}_metabin2/
cp 05.binning/results/02.concoct.metahit.bins/${i}.fa 05.binning/bins_metahit/${i}_concoct/
Fasta_to_Scaffolds2Bin.sh -i 05.binning/bins_metahit/${i}_metabat2/ -e fa > 06.das_tool/${i}_metahit_metabat2.scaffolds2bin.tsv
Fasta_to_Scaffolds2Bin.sh -i 05.binning/bins_metahit/${i}_metabin2/ -e fasta > 06.das_tool/${i}_metahit_metabin2.scaffolds2bin.tsv
Fasta_to_Scaffolds2Bin.sh -i 05.binning/bins_metahit/${i}_concoct/ -e fa > 06.das_tool/${i}_metahit_concoct.scaffolds2bin.tsv
done

for i in $(cat data/all_sampe_id); do
cp 05.binning/results/00.metabat2.metaspades.bins/${i}*.fa 05.binning/bins_metaspades/${i}_metabat2/
cp 05.binning/results/01.metabin2.metaspades.bins/${i}_metabin2.*.fasta 05.binning/bins_metaspades/${i}_metabin2/
cp 05.binning/results/02.concoct.metaspades.bins/${i}.fa 05.binning/bins_metaspades/${i}_concoct/
Fasta_to_Scaffolds2Bin.sh -i 05.binning/bins_metaspades/${i}_metabat2/ -e fa > 06.das_tool/${i}_metaspades_metabat2.scaffolds2bin.tsv
Fasta_to_Scaffolds2Bin.sh -i 05.binning/bins_metaspades/${i}_metabin2/ -e fasta > 06.das_tool/${i}_metaspades_metabin2.scaffolds2bin.tsv
Fasta_to_Scaffolds2Bin.sh -i 05.binning/bins_metaspades/${i}_concoct/ -e fa > 06.das_tool/${i}_metaspades_concoct.scaffolds2bin.tsv
done
# dereplication
##bins generated based on metahit assembly
for i in $(cat data/all_sampe_id); do
DAS_Tool -i 06.das_tool/${i}_metahit_metabat2.scaffolds2bin.tsv,06.das_tool/${i}_metahit_metabin2.scaffolds2bin.tsv ${i}_metahit_concoct.scaffolds2bin.tsv, -l metabat,maxbin,concoct -c 04.assmebly_firlter/00.assembly_filter_results_megahit/${i}.fa -o 07.das_tool_results/metahit/${i} --search_engine diamond --write_bins 1 --score_threshold 0 -t 30
done
##bins generated based on metaspades assembly
for i in $(cat data/all_sampe_id); do
DAS_Tool -i 06.das_tool/${i}_metaspades_metabat2.scaffolds2bin.tsv,06.das_tool/${i}_metaspades_metabin2.scaffolds2bin.tsv 06.das_tool/${i}_metaspades_concoct.scaffolds2bin.tsv, -l metabat,maxbin,concoct -c 04.assmebly_firlter/01.assembly_filter_results_metspades/${i}.fa -o 07.das_tool_results/metaspades/${i} --search_engine diamond --write_bins 1 --score_threshold 0 -t 30
done
```

## 6.CheckM

```sh
for i in $(cat data/all_sampe_id); do
checkm lineage_wf -x fa 07.das_tool_results/metahit/${i}_DASTool_bins/ 08.checkm/${i}.metahit -t 30 --tmpdir bin_checkm.tmp
checkm lineage_wf -x fa 07.das_tool_results/metaspades/${i}_DASTool_bins/ 08.checkm/${i}.metaspades -t 30 --tmpdir bin_checkm.tmp
done
```

## 7.dRep

```sh
## 99% ANI (strain level)
for i in Cheese Chinese_liquor Cocoa Coffee Ganjang Kimchi Kombucha Koumiss Lajiaojiang Milk_kefir Nunu Sour_dough Thick_broad_bean_sauce Wine Yoghurt; do
dRep dereplicate 10.dRep/drep_${i}/ -g 09.genome/${i}/*.fa -p 2 --ignoreGenomeQuality -pa 0.95 -sa 0.99 --S_algorithm fastANI
done
## 95% ANI (species level)
for i in Cheese Chinese_liquor Cocoa Coffee Ganjang Kimchi Kombucha Koumiss Lajiaojiang Milk_kefir Nunu Sour_dough Thick_broad_bean_sauce Wine Yoghurt; do
dRep dereplicate -g 10.dRep/drep_${i}/dereplicated_genomes/*.fa -p 20 11.SGB/${i}/ --ignoreGenomeQuality -pa 0.9 -sa 0.95 --S_algorithm fastANI
done
```

## 8. GTDB_tk

```sh
#Developing a MAGs dataset
for i in Cheese Chinese_liquor Cocoa Coffee Ganjang Kimchi Kombucha Koumiss Lajiaojiang Milk_kefir Nunu Sour_dough Thick_broad_bean_sauce Wine Yoghurt; do
cp 10.dRep/drep_${i}/dereplicated_genomes/*.fa 12.GTDB_tk/00.genome/
done
##Taxon classification
gtdbtk classify_wf --genome_dir 12.GTDB_tk/00.genome/ --out_dir 12.GTDB_tk/01.genome_gtdbtk --cpus 30 --scratch_dir gtdbtk.temp --pplacer_cpus 32 --extension fa
```

## 9. Phylogenetic analysis

```sh
##Phylogenetic analysis
bmge -i 12.GTDB_tk/01.genome_gtdbtk/gtdbtk.bac120.user_msa.fasta -t AA -g 0.5 -h 1 -b 1 -w 1 -of 13.tree/out_bmge_trimmed_2.fasta
FastTree out_bmge_trimmed_2.fasta > 653MAGs.tree
```

## 10.  antismash analysis

```sh
for i in $(cat data/MAG_id); do
antismash 12.GTDB_tk/00.genome/${i}.fa --taxon bacteria --output-dir 13.SM_BGC/01.antismash_fullhmmer/${i} --genefinding-tool prodigal --cb-knownclusters -c 3 --cc-mibig --fullhmmer
done
```

## 11. Bio-SCAPE

```sh
##Preparing configuration Files
for i in $(cat data/MAG_id); do
cp 13.SM_BGC/01.antismash_fullhmmer/${i}/*region*gbk 13.SM_BGC/input_gbk
done

python ./BiG-SCAPE-master/bigscape.py -i 13.SM_BGC/input_gbk -o output_bigscape --pfam_dir ./BiG-SCAPE-master/hmm -c 30 --cutoffs 0.7 --include_singletons --mode auto --mix --no_classify --mibig --verbose
```

## 12. Bio-SLiCE

```sh
##BiG-SLiCE clustering analysis, the input folder was prepared as previously reported (https://github.com/medema-group/bigslice)
bigslice -i INPUT OUTPUT

##-query mode
bigslice --query 13.SM_BGC/input_gbk --n_ranks 1 OUTPUT -t 30
```

## 13. BGC_function_prediction

```sh
for i in $(cat data/BGC_id); do
./cluster_function_prediction.py 13.SM_BGC/input_gbk/${i}.gbk 14.BGC_function/rgi/${i}.txt --output 14.BGC_function/BGC_function_prediction --antismash_version 5 --rgi_version 5
done
```

