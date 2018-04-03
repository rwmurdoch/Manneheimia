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



