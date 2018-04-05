# Manneheimia
Description of pipeline for processing of Mannheimia amplicon library

# Overview

This is a description of a pipeline that was used for processing a set of ~180 amplicon libraries.  

## 1. qiime2_master_RWM.sh 

The user first must download and install the qiime2 2018.2 package (you're on your own with this).

Open and read through the script.  There are a number of things that must be altered for any given project, including project name, data locations, reference database location and indentity, primer sequences, sample filtering thresholds, and CPU thread counts.  All important information along these lines is explained on # comment lines in the .sh file.

This script contains all of the qiime2 code needed to:
1. Import and configure a Silva reference database
2. Import a paired-end demultiplexed set of amplicon libraries
3. Import a metadata file
4. Run the project through three different pipelines, all of which will create an OTU table, perform taxonomic placement, and generate stacked bar charts.  It is assumed that the user will check the results at this stage and choose a pipeline to move forward with.
5. Filter samples based on a % of median approach, which the user must generate indepently. This can be done by getting seq counts from the appropriate qzv files.
6. Generate full diversity metrics, alpha and beta diversity, ordinations, and corresponding qzv visualization files

At this stage, the user can visualize any .qzv file at the qiime2view website.  This pipeline takes the:
1. OTU table
2. taxonomy table
3. rooted tree

and transitions into the phyloseq R package for further analysis.

### explicit pipeline details:

cutadapt was used to detect and remove primer sequences using default parameters
Taxonomic placement: Naive-Bayes feature classifier (sklearn) trained on Silva 128 reference database filtered using the primers utilized in the study.

#### VSearch open reference clustering

1. join pairs, 50 bp minimum overlap, 380 minimum and 480 maximum merged sequence length, disposing of sequences that do not meet these specifications.
2. quality score read filtering, q score 25 minimum with default settings
3. de novo chimera searching using uchime
4. vsearch open reference clustering using silva v.128 97% database
5. taxonomic placement using trained sklearn model

#### deblur16S

1. join pairs, 50 bp minimum overlap, 380 minimum and 480 maximum merged sequence length, disposing of sequences that do not meet these specifications.
2. deblur denoise-16S for exact sequence variant identification and filtering of non-16S sequences
3. taxonomic placement using trained sklearn model

#### dada2

1. raw data fed directly to data2 denoise-paired module for exact sequence variants, truncating forward read to 250bp and reverse read to 200bp
2. taxonomic placement using trained sklearn model

#### further processing of the dada2 results

1. remove all samples with read counts below 5% of the median (1200)
2. align using MAFFT
3. midpoint-rooted tree using fasttree


## 2. preparing qiime2 files for importing to physeq

moving the qiime2 results into physeq is a bit laborious.  Qiime2 is new and the hand-off isn't smooth yet.  Ways to handle this are discussed here: https://github.com/joey711/phyloseq/issues/821 

The explicit instructions for getting data out of qiime2 .qza files and into physeq begin below (taken from the link above):

based on this thread: https://forum.qiime2.org/t/converting-biom-files-with-taxonomic- info-for-import-in-r-with-phyloseq/2542/5

### How to export a feature (OTU) table and convert from biom to .tsv (use following codes in QIIME2)
#### Step 1, export OTU table
qiime tools export \ yourOTUtablename.qza \ --output-dir phyloseq

#### OTU table exports as feature-table.biom so convert to .tsv
biom convert -i phyloseq/feature-table.biom -o phyloseq/otu_table.txt --to-tsv 

now you have otu_table.txt
open it up in text edit and change #OTUID to OTUID
#### Step 2, export taxonomy table
qiime tools export \ taxonomy.qza\
--output-dir phyloseq
now you have taxonomy.tsv
open it up in text edit and change Feature ID to OTUID
#### Step 3, export tree
qiime tools export \ unrooted-tree.qza \ --output-dir phyloseq
#### Step 4, if you filtered out any sequences (chloroplasts, mitochondria, etc) then your taxonomy and OTU tables are different lengths. QIIME2 doesn’t filter out taxonomy, so you have to merge the two files in R and output a merged file.
(use following codes in R)

### qiime2_to_physeq1.R

This first script will take the qiime2 files and begin to format them for physeq.  After adapting it to your project, you will have to do some work by hand:

### manually using Excel to format the tables

It seems tedious but you need to open the merged .txt file in excel and split into two files: one for taxonomy (containing only the columns OTUID and taxonomic info) and the other for the OTU matrix (containing only OTUID and abundances in each sample). Note: for the taxonomy file, you need to use data —> text-to-columns in Excel and separate on semicolon to get columns for kingdom, phylum, class, etc... once you make these two separate files in excel, save each as a .csv

### qiime2_to_physeq.R

This second R script will handle taking the files produced above and convert them into physeq objects (within R)

## split_and_properties.R

In this project, we will first make sure that the sub-5% threshold samples have been removed.  Then, the script splits the project by sampling location.

http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531
https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis

The two references above address an important primary question; should you normalize or rarify your data?  The best answer seems to be no, no normalization is appropriate.  However, sample depth should be considered.

1. Split by samplling location
2. generate read-depth figures
3. calculate alpha diversity plots
