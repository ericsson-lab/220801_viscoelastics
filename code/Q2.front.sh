#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH --job-name=Q2.front.sh       # name for the job
#SBATCH --cpus-per-task=3              # number of cores
#SBATCH --mem=100G                       # total memory
#SBATCH --nodes 2
#SBATCH --time 2-00:00:00                 # time limit in the form days-hours:minutes
#SBATCH --mail-user=zlmg2b@umsystem.edu    # email address for notifications
#SBATCH --mail-type=FAIL,END,BEGIN           # email types
#SBATCH --partition Lewis            # max of 1 node and 4 hours; use `Lewis` for larger jobs
#--------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

module load miniconda3
source activate qiime2-2021.8

mkdir sequences
mkdir $4
cd $4
mkdir Dada2
mkdir Dada2/visualization
mkdir phylogeny
mkdir sequences
mkdir taxonomy
mkdir transfer
mkdir transfer/sequences
mkdir transfer/Dada2
mkdir transfer/taxonomy
cd ..

# python3 ../workflow/generate_manifest.py -i $1 -l $2 -m $3

cd demux_seqs

for file in *R1*.gz; do [ -f "$file" ] || continue; mv -vf "$file" "${file//_*R1_001.fastq.gz/_R1.fastq.gz}"; done
for file in *R2*.gz; do [ -f "$file" ] || continue; mv -vf "$file" "${file//_*R2_001.fastq.gz/_R2.fastq.gz}"; done

cd ..

cp ./metadata.txt ./$4/metadata.txt

## Import demuxed seq (*.fastq.gz) from manifest file
## manifest generated in generate_manifest.py

qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33 \
  --input-path ./manifest.csv \
  --output-path ./sequences/demux_seqs.qza

cp ./manifest.csv ./$4/manifest.csv

## Visualize imported seqs
qiime demux summarize \
  --i-data ./sequences/demux_seqs.qza \
  --o-visualization ./sequences/demux_seqs.qzv

qiime tools export \
  --input-path ./sequences/demux_seqs.qzv \
  --output-path ./sequences/

mv ./sequences/per-sample-fastq-counts.tsv \
   ./sequences/per-sample-fastq-counts_untrimmed.tsv

cd ./sequences/
rm -rf data.jsonp
rm -rf demultiplex-summary*
rm -rf dist
rm -rf *.html
rm -rf q2templateassets
rm -rf *seven-number-summaries.tsv
cd ..

## Trim Ilumina primers/adapters from demux seqs
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ./sequences/demux_seqs.qza \
  --p-cores $SLURM_CPUS_ON_NODE \
  --p-adapter-f 'ATTAGAWACCCBDGTAGTCC' \
  --p-front-f 'GTGCCAGCMGCCGCGGTAA' \
  --p-adapter-r 'TTACCGCGGCKGCTGGCAC' \
  --p-front-r 'GGACTACHVGGGTWTCTAAT' \
  --p-discard-untrimmed \
  --p-no-indels \
  --o-trimmed-sequences ./sequences/trimmed_demux_seqs.qza 

## Visualize trimmed seqs
qiime demux summarize \
  --i-data ./sequences/trimmed_demux_seqs.qza \
  --o-visualization ./sequences/trimmed_demux_seqs.qzv

qiime tools export \
  --input-path ./sequences/trimmed_demux_seqs.qzv \
  --output-path ./sequences/

mv ./sequences/per-sample-fastq-counts.tsv \
   ./sequences/per-sample-fastq-counts_trimmed.tsv

cd ./sequences/
rm -rf data.jsonp
rm -rf demultiplex-summary*
rm -rf dist
rm -rf *.html
rm -rf q2templateassets
rm -rf *seven-number-summaries.tsv
cd ..

cp ./sequences/*.tsv ./$4/sequences/

cd $4

## Denoise seqs to unique Amplicon Sequence Variants (ASVs)
## Parameters are historical values from IRCF pipeline
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ../sequences/trimmed_demux_seqs.qza \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --o-table ./Dada2/dada2_table.qza \
  --o-representative-sequences ./Dada2/dada2_rep_seqs.qza \
  --o-denoising-stats ./Dada2/dada2_stats.qza \
  --p-n-threads $SLURM_CPUS_ON_NODE

## Visualize Dada2 stats to assess filtering results
qiime metadata tabulate \
  --m-input-file ./Dada2/dada2_stats.qza  \
  --o-visualization ./Dada2/visualization/dada2_stats.qzv

## Filter seqs down to 249-257 in length (based on historical pipeline)
qiime feature-table filter-seqs \
  --i-data ./Dada2/dada2_rep_seqs.qza \
  --m-metadata-file ./Dada2/dada2_rep_seqs.qza \
  --p-where 'length(sequence) >= 249 AND length(sequence) <= 257' \
  --o-filtered-data ./Dada2/dada2_rep_seqs_filtered.qza 

## Filter feature table based on retained seqs
qiime feature-table filter-features \
  --i-table ./Dada2/dada2_table.qza \
  --m-metadata-file ./Dada2/dada2_rep_seqs_filtered.qza \
  --o-filtered-table ./Dada2/dada2_table_filtered.qza

## Visualize feature table using seqs filtered to 249-257 bp
qiime feature-table summarize \
  --i-table ./Dada2/dada2_table_filtered.qza \
  --m-sample-metadata-file ../metadata.txt \
  --o-visualization ./Dada2/visualization/dada2_table_filtered.qzv

## Export Dada2 counts per sample per feature CSV files 
qiime tools export \
  --input-path ./Dada2/visualization/dada2_table_filtered.qzv \
  --output-path ./Dada2/visualization/

qiime tools export \
  --input-path ./Dada2/visualization/dada2_stats.qzv \
  --output-path ./Dada2/visualization/

## Cleanup Dada2 export 
rm -rf ./Dada2/visualization/*.html
rm -rf ./Dada2/visualization/*.pdf
rm -rf ./Dada2/visualization/*.png
rm -rf ./Dada2/visualization/js
rm -rf ./Dada2/visualization/licenses
rm -rf ./Dada2/visualization/q2templateassets
rm -rf ./Dada2/visualization/css

## rename Dada2 feature table csv files
mv ./Dada2/visualization/feature-frequency-detail.csv  \
   ./Dada2/visualization/dada2-feature-frequency-detail.csv
mv ./Dada2/visualization/sample-frequency-detail.csv \
   ./Dada2/visualization/dada2-sample-frequency-detail.csv
#mv ./Dada2/visualization/sample-frequency-detail.csv \
#   ./Dada2/visualization/dada2-sample-frequency-detail.csv

## Phylogeny generated using denovo approach
## Trees used for phylogenetic diversity metrics and empress plots 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ./Dada2/dada2_rep_seqs_filtered.qza \
  --o-alignment ./phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment ./phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree ./phylogeny/unrooted-tree.qza \
  --o-rooted-tree ./phylogeny/rooted-tree.qza

## Assign taxonomy based on Silva 138 ref database
qiime feature-classifier classify-sklearn \
  --i-reads ./Dada2/dada2_rep_seqs_filtered.qza \
  --i-classifier /home/zlmg2b/data/workflow/feature-classifier/silva-138-99-classifier-515-806.qza \
  --o-classification ./taxonomy/taxonomy.qza

## Export taxonomy to TSV file
qiime tools export \
  --input-path ./taxonomy/taxonomy.qza \
  --output-path ./taxonomy/

## Matches taxonomy to sequence
qiime metadata tabulate \
  --m-input-file ./Dada2/dada2_rep_seqs_filtered.qza \
  --m-input-file ./taxonomy/taxonomy.qza \
  --o-visualization ./Dada2/visualization/dada2_rep_seqs_filtered_taxa.qzv

## Export dada2_rep_seqs_filtered_taxa.qzv for results file
qiime tools export \
  --input-path ./Dada2/visualization/dada2_rep_seqs_filtered_taxa.qzv \
  --output-path ./taxonomy/
## Rename file
mv ./taxonomy/metadata.tsv ./taxonomy/dada2_rep_seqs_filtered_taxa.tsv

## Remove excess files from exporitng dada2_rep_seqs_filtered_taxa.qzv
rm -r ./taxonomy/index.html
rm -r ./taxonomy/q2templateassets
rm -r ./taxonomy/js
rm -r ./taxonomy/css

## Exports Dada2 feature table and renames it
qiime tools export \
  --input-path ./Dada2/dada2_table_filtered.qza \
  --output-path ./Dada2/
mv ./Dada2/feature-table.biom ./Dada2/dada2_feature-table.biom

## Adds taxonomy labels to Dada2 feature table
biom add-metadata \
  -i Dada2/dada2_feature-table.biom \
  -o Dada2/dada2_feature-table_taxa.biom \
  --observation-metadata-fp taxonomy/taxonomy.tsv \
  --observation-header="Feature ID,Taxon" \
  --sc-separated taxonomy

## Converts to CSV file
biom convert \
  -i Dada2/dada2_feature-table_taxa.biom \
  -o Dada2/dada2_feature-table_taxa.tsv \
  --to-tsv \
  --output-metadata-id=Taxon \
  --tsv-metadata-formatter=naive \
  --header-key=Taxon

cp ../sequences/per-sample-fastq-counts_untrimmed.tsv \
   ./transfer/sequences/per-sample-fastq-counts_untrimmed.tsv
cp ../sequences/per-sample-fastq-counts_untrimmed.tsv \
   ./transfer/sequences/per-sample-fastq-counts_trimmed.tsv
cp ../sequences/demux_seqs.qzv ./transfer/sequences/demux_seqs.qzv
cp ../sequences/trimmed_demux_seqs.qzv ./transfer/sequences/trimmed_demux_seqs.qzv
cp ./Dada2/visualization/dada2_stats.qzv ./transfer/Dada2/dada2_stats.qzv
cp ./Dada2/visualization/dada2_table_filtered.qzv ./transfer/Dada2/dada2_table_filtered.qzv
cp ./taxonomy/taxonomy.tsv ./transfer/taxonomy/taxonomy.tsv
cp ./Dada2/visualization/dada2_rep_seqs_filtered_taxa.qzv \
   ./transfer/taxonomy/dada2_rep_seqs_filtered_taxa.qzv
cp ./Dada2/dada2_feature-table_taxa.tsv ./transfer/Dada2/dada2_feature-table_taxa.tsv
cp ./taxonomy/dada2_rep_seqs_filtered_taxa.tsv ./transfer/taxonomy/dada2_rep_seqs_filtered_taxa.tsv
cp ./Dada2/visualization/dada2-sample-frequency-detail.csv \
   ./transfer/Dada2/dada2-sample-frequency-detail.csv
cp ./Dada2/visualization/dada2-feature-frequency-detail.csv \
   ./transfer/Dada2/dada2-feature-frequency-detail.csv
cp ./taxonomy/taxonomy.qza \
   ./transfer/taxonomy/taxonomy.qza
mv ./metadata.txt ./metadata_$4.txt
cp ./phylogeny/rooted-tree.qza \
   ./transfer/

## TO ADD
## Export Dada2 stats
## Export read counts
## generate Excel file
## remove line 1 of tables
## change #OTUID to featureid in those tables


echo "### Ending at: $(date) ###"