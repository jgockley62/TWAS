---
title: Covariate and analysis of MSSM reprocessed counts Regressing out Diagnosis as well for TWAS
author: 'JKG - Adapted from: Thanneer Perumal'
output: html_notebook
editor_options:
  chunk_output_type: console
---
  Date of analysis update: "Thu Sep 19 00:32:33 2019"
  
  
syn8691099.1
syn10156693.3
syn8698270.1
syn8449369.2
syn6100548.7
syn6101474.4





```
## Welcome, Jake Gockley!
```

```
## NULL
```

### Data download
Obtain count matrix and metadata from synapse

```
## Downloading  [#-------------------]4.28%   2.0MB/46.7MB (54.2MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [##------------------]8.56%   4.0MB/46.7MB (64.8MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [###-----------------]12.84%   6.0MB/46.7MB (73.0MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [###-----------------]17.12%   8.0MB/46.7MB (81.0MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [####----------------]21.40%   10.0MB/46.7MB (87.8MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [#####---------------]25.69%   12.0MB/46.7MB (92.9MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [######--------------]29.97%   14.0MB/46.7MB (96.9MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [#######-------------]34.25%   16.0MB/46.7MB (99.9MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [########------------]38.53%   18.0MB/46.7MB (102.5MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [#########-----------]42.81%   20.0MB/46.7MB (104.9MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [#########-----------]47.09%   22.0MB/46.7MB (106.8MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [##########----------]51.37%   24.0MB/46.7MB (109.2MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [###########---------]55.65%   26.0MB/46.7MB (110.6MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [############--------]59.93%   28.0MB/46.7MB (111.3MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [#############-------]64.21%   30.0MB/46.7MB (112.3MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [##############------]68.50%   32.0MB/46.7MB (114.5MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [###############-----]72.78%   34.0MB/46.7MB (116.1MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [###############-----]77.06%   36.0MB/46.7MB (117.8MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [################----]81.34%   38.0MB/46.7MB (119.0MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [#################---]85.62%   40.0MB/46.7MB (119.9MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [##################--]89.90%   42.0MB/46.7MB (121.1MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [###################-]94.18%   44.0MB/46.7MB (122.3MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [####################]98.46%   46.0MB/46.7MB (123.8MB/s) MSSM_all_counts_matrix.txt.gz     Downloading  [####################]100.00%   46.7MB/46.7MB (124.7MB/s) MSSM_all_counts_matrix.txt.gz Done...
```

```
## Downloading  [####################]100.00%   294.9kB/294.9kB (27.8MB/s) MSBB_RNAseq_covariates.csv Done...
```

```
## Downloading  [####################]100.00%   8.9kB/8.9kB (16.2MB/s) MSBB_Age90andAbove.tsv Done...
```

```
## Downloading  [####################]100.00%   20.8kB/20.8kB (31.0MB/s) MSBB_clinical.csv Done...
```

```
## Downloading  [####################]100.00%   373.8kB/373.8kB (55.4MB/s) MSSM_all_metrics_matrix.txt Done...
```

```
## Downloading  [#############-------]67.05%   2.0MB/3.0MB (75.1MB/s) geneParameters.tsv     Downloading  [####################]100.00%   3.0MB/3.0MB (39.6MB/s) geneParameters.tsv Done...
```

### Data preprocessing

754 samples from 253 subjects were obtained from the MSSM cohorts in AMP-AD reprocessed rnaseq project. 

Following 272 samples were removed: BM_10_569, BM_10_571, BM_10_573, BM_10_594, BM_10_621, BM_10_634, BM_10_652, BM_10_677, BM_10_741, BM_10_750, BM_10_751, BM_10_754, BM_10_766, BM_10_798, BM_22_125, BM_22_131, BM_22_132, BM_22_134, BM_22_135, BM_22_150, BM_22_160, BM_22_211, BM_22_258, BM_22_263, BM_22_41, BM_22_44, BM_22_49, BM_22_57, BM_22_64, BM_36_308, BM_36_327, BM_36_338, BM_36_344, BM_36_347, BM_36_348, BM_36_349, BM_36_390, BM_36_429, BM_36_432, BM_36_485, BM_36_491, BM_36_492, BM_36_495, BM_36_505, BM_36_516, BM_36_527, BM_36_533, BM_36_539, BM_36_541, hB_RNA_10232, hB_RNA_10242, hB_RNA_10252, hB_RNA_10362, hB_RNA_10402, hB_RNA_10412, hB_RNA_10482_resequenced, hB_RNA_10492_resequenced, hB_RNA_10502, hB_RNA_10502_resequenced, hB_RNA_10522, hB_RNA_10522_resequenced, hB_RNA_10532, hB_RNA_10532_resequenced, hB_RNA_10542, hB_RNA_10542_resequenced, hB_RNA_10552, hB_RNA_10567, hB_RNA_10567_resequenced, hB_RNA_10577, hB_RNA_10577_resequenced, hB_RNA_10583, hB_RNA_10617, hB_RNA_10632, hB_RNA_10642, hB_RNA_10652, hB_RNA_10652_resequenced, hB_RNA_10662, hB_RNA_10672, hB_RNA_10682, hB_RNA_10692, hB_RNA_10692_resequenced, hB_RNA_10712, hB_RNA_10712_resequenced, hB_RNA_10722_resequenced, hB_RNA_10742_resequenced, hB_RNA_10762_resequenced, hB_RNA_10782, hB_RNA_10802, hB_RNA_10802_resequenced, hB_RNA_10822, hB_RNA_10832, hB_RNA_10842, hB_RNA_10852, hB_RNA_10852_resequenced, hB_RNA_10862, hB_RNA_10862_resequenced, hB_RNA_10872, hB_RNA_10872_resequenced, hB_RNA_10882, hB_RNA_10902, hB_RNA_10912, hB_RNA_10922, hB_RNA_10932, hB_RNA_10942, hB_RNA_10952, hB_RNA_11023, hB_RNA_11082, hB_RNA_11092, hB_RNA_11102, hB_RNA_11112, hB_RNA_11122, hB_RNA_11152, hB_RNA_12171, hB_RNA_12171_resequenced, hB_RNA_12181, hB_RNA_12181_resequenced, hB_RNA_12191, hB_RNA_12201, hB_RNA_12211, hB_RNA_12222, hB_RNA_12232, hB_RNA_12242, hB_RNA_12252, hB_RNA_12252_resequenced, hB_RNA_12262, hB_RNA_12262_resequenced, hB_RNA_12272, hB_RNA_12272_resequenced, hB_RNA_12282, hB_RNA_12292, hB_RNA_12312, hB_RNA_12312_resequenced, hB_RNA_12322, hB_RNA_12342, hB_RNA_12342_resequenced, hB_RNA_12352, hB_RNA_12352_resequenced, hB_RNA_12372, hB_RNA_12372_resequenced, hB_RNA_12382, hB_RNA_12382_resequenced, hB_RNA_12402, hB_RNA_12464, hB_RNA_12484, hB_RNA_12514, hB_RNA_12532, hB_RNA_12555, hB_RNA_12606, hB_RNA_12615, hB_RNA_12704, hB_RNA_12734, hB_RNA_12762, hB_RNA_12784, hB_RNA_12838, hB_RNA_12859, hB_RNA_12887, hB_RNA_12964, hB_RNA_13032, hB_RNA_13048, hB_RNA_13048_resequenced, hB_RNA_13058, hB_RNA_13058_resequenced, hB_RNA_13068_resequenced, hB_RNA_13081, hB_RNA_13081_resequenced, hB_RNA_13107, hB_RNA_13134, hB_RNA_13144, hB_RNA_13216_resequenced, hB_RNA_13228, hB_RNA_13253, hB_RNA_13266, hB_RNA_13289, hB_RNA_13359, hB_RNA_13389, hB_RNA_13518, hB_RNA_13547, hB_RNA_13649, hB_RNA_16295, hB_RNA_16495, hB_RNA_16635, hB_RNA_16805, hB_RNA_16815, hB_RNA_17045, hB_RNA_4301, hB_RNA_4341, hB_RNA_4411, hB_RNA_4429, hB_RNA_4501, hB_RNA_4531, hB_RNA_4546, hB_RNA_4591, hB_RNA_4697, hB_RNA_4751, hB_RNA_4791, hB_RNA_4801, hB_RNA_4923, hB_RNA_4946, hB_RNA_4981, hB_RNA_4991, hB_RNA_5031, hB_RNA_5161, hB_RNA_7735, hB_RNA_7765_resequenced, hB_RNA_7785, hB_RNA_7805, hB_RNA_7825, hB_RNA_7855, hB_RNA_8075_resequenced, hB_RNA_8085_resequenced, hB_RNA_8105, hB_RNA_8115_resequenced, hB_RNA_8125_resequenced, hB_RNA_8155_resequenced, hB_RNA_8175, hB_RNA_8175_resequenced, hB_RNA_8215_resequenced, hB_RNA_8255, hB_RNA_8265, hB_RNA_8265_resequenced, hB_RNA_8285, hB_RNA_8285_resequenced, hB_RNA_8295, hB_RNA_8315, hB_RNA_8335, hB_RNA_8345, hB_RNA_8365, hB_RNA_8395, hB_RNA_8405, hB_RNA_8425, hB_RNA_8445, hB_RNA_8475, hB_RNA_8575, hB_RNA_8645, hB_RNA_8685, hB_RNA_8705, hB_RNA_8795, hB_RNA_8835, hB_RNA_8915, hB_RNA_9045, hB_RNA_9072, hB_RNA_9074, hB_RNA_9085, hB_RNA_9105, hB_RNA_9115, hB_RNA_9125, hB_RNA_9137, hB_RNA_9141, hB_RNA_9151, hB_RNA_9158, hB_RNA_9162, hB_RNA_9168, hB_RNA_9170, hB_RNA_9171, hB_RNA_9173, hB_RNA_9174, hB_RNA_9176, hB_RNA_9191, hB_RNA_9194, hB_RNA_9207, hB_RNA_9207_resequenced, hB_RNA_9208, hB_RNA_9208_resequenced, hB_RNA_9209_resequenced, hB_RNA_9210_resequenced, hB_RNA_9212_resequenced, hB_RNA_9213, hB_RNA_9215, hB_RNA_9216, hB_RNA_9222, hB_RNA_9223, hB_RNA_9225

### Covariate clustering
Determine relationship between covariates

Correlation/association between covariates at an FDR <= 0.1
${image?fileName=covariates%2Ecorrelation%2D1%2Epng&align=none&scale=100}
### Explore metadata
${image?fileName=data%2Eexplore%2D1%2Epng&align=none&scale=100}

### Filter genes
* Remove genes that have less than 1 cpm counts in at least 50% of samples per Tissue x Diagnosis
* Remove genes with missing gene length and percentage GC content

|Biotype                            |  fraction|
|:----------------------------------|---------:|
|antisense                          | 0.0438028|
|lincRNA                            | 0.0347486|
|processed_pseudogene               | 0.0159060|
|protein_coding                     | 0.8468739|
|transcribed_unprocessed_pseudogene | 0.0118683|

### Library Normalisation
Library normalisation is performed using cqn (conditional quantile normalisation)

```
## fitting ...
##   |                                                                         |                                                                 |   0%  |                                                                         |======================                                           |  33%  |                                                                         |===========================================                      |  67%  |                                                                         |=================================================================| 100%
```

### Outlier Analysis
#### Sample outliers
Outlier analysis is performed before library normalisation with raw cpm counts
${image?fileName=outlier%2Eanalysis%2D1%2Epng&align=none&scale=100}
${image?fileName=outlier%2Eanalysis%2D2%2Epng&align=none&scale=100}
Processing 16346 genes in 753 samples

Based on the expression pattern following samples were tagged as outliers: hB_RNA_10622

Distribution of samples are: 

|Tissue | AD| CONTROL| OTHER|
|:------|--:|-------:|-----:|
|FP     | 90|      45|    79|
|IFG    | 79|      37|    71|
|PHG    | 65|      38|    58|
|STG    | 85|      37|    69|



#### Gene outliers
Assign NA values to genes that are above and below 3 std deviation of its distribution

### Sample clustering
PCA based clustering of samples
${image?fileName=decompse%2Enormalise%2Edata1%2E1%2D1%2Epng&align=none&scale=100}
Tree based classification of samples
${image?fileName=decompse%2Enormalise%2Edata1%2E2%2D1%2Epng&align=none&scale=100}
Coexpression of genes 
${image?fileName=coexp1%2D1%2Epng&align=none&scale=100}

### Significant Covariates
Correlation between pca of unadjusted mRNA expression and covariates are used to find significant covariates

```
## 
## Running PCA and calculating correlations for:
## Scaled NULL design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations
```

Significant covariates to adjust at FDR 0.1 are individualIdentifier, batch, RACE, Sex, Tissue.Diagnosis, Tissue.APOE4, Diagnosis, CDR, RIN, RIN2, AOD, PMI, PCT_PF_READS_ALIGNED, PCT_CODING_BASES, PCT_INTERGENIC_BASES, PCT_INTRONIC_BASES, PCT_RIBOSOMAL_BASES
${image?fileName=preadj%2Ecovariates%2Eplot%2D1%2Epng&align=none&scale=100}

### Normalisation (iterative design)
Since many covariates are correlated, re-normalising and re-adjusting COUNTS with an iterative design matrix.

NOTE:
1. Using a mixed effect model where random effect is chosen as individualIdentifier
2. Adding batch and Sex a priori to variable selection
3. Primary variable of interest Tissue.Diagnosis is excluded from the pool of available covariates for selection
Using following covariates in the model: batch, Sex, PCT_INTRONIC_BASES as fixed effects and individualIdentifier as random effect

Running PCA and calculating correlations for:
Scaled adjusted design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations 
Using following covariates in the model: batch, Sex, PCT_INTRONIC_BASES, RIN as fixed effects and individualIdentifier as random effect

Running PCA and calculating correlations for:
Scaled adjusted design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations 
Using following covariates in the model: batch, Sex, PCT_INTRONIC_BASES, RIN, PMI as fixed effects and individualIdentifier as random effect

Running PCA and calculating correlations for:
Scaled adjusted design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations 
Using following covariates in the model: batch, Sex, PCT_INTRONIC_BASES, RIN, PMI, PCT_INTERGENIC_BASES as fixed effects and individualIdentifier as random effect

Running PCA and calculating correlations for:
Scaled adjusted design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations 
Using following covariates in the model: batch, Sex, PCT_INTRONIC_BASES, RIN, PMI, PCT_INTERGENIC_BASES, AOD as fixed effects and individualIdentifier as random effect

Running PCA and calculating correlations for:
Scaled adjusted design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations 
Using following covariates in the model: batch, Sex, PCT_INTRONIC_BASES, RIN, PMI, PCT_INTERGENIC_BASES, AOD, Diagnosis as fixed effects and individualIdentifier as random effect

Running PCA and calculating correlations for:
Scaled adjusted design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations 
Using following covariates in the model: batch, Sex, PCT_INTRONIC_BASES, RIN, PMI, PCT_INTERGENIC_BASES, AOD, Diagnosis, PCT_CODING_BASES as fixed effects and individualIdentifier as random effect

Running PCA and calculating correlations for:
Scaled adjusted design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations 
Using following covariates in the model: batch, Sex, PCT_INTRONIC_BASES, RIN, PMI, PCT_INTERGENIC_BASES, AOD, Diagnosis, PCT_CODING_BASES, RACE as fixed effects and individualIdentifier as random effect

Running PCA and calculating correlations for:
Scaled adjusted design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations 
Using following covariates in the final model: batch, Sex, PCT\_INTRONIC\_BASES, RIN, PMI, PCT\_INTERGENIC\_BASES, AOD, Diagnosis, PCT\_CODING\_BASES, RACE, PCT\_PF\_READS\_ALIGNED as fixed effects

### Sanity check

```
## 
## Running PCA and calculating correlations for:
## Scaled adjusted design(voom-normalized)  data in PCA; PVE >= 1%; pearson correlations
```

${image?fileName=residual%2Eadj%2D1%2Epng&align=none&scale=100}
Coexpression of genes 
${image?fileName=coexp2%2D1%2Epng&align=none&scale=100}
PCA of residual data
${image?fileName=decompse%2Enormalise%2Edata2%2E1%2D1%2Epng&align=none&scale=100}
Tree based clustering of residual data
${image?fileName=decompse%2Enormalise%2Edata2%2E2%2D1%2Epng&align=none&scale=100}

### Adjust data with covariates for Network Analysis
Identified covariates are regressed out from the expression matrix for network analysis


### SVA Adjustments for eQTL analysis
Conditioned on primary variable (Tissue.Diagnosis) and identified covariates estimate surrogate variables using SVA package
${image?fileName=sva%2Eadjust%2D1%2Epng&align=none&scale=100}

```
## Number of significant surrogate variables is:  18 
## Iteration (out of 30 ):1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30
```






