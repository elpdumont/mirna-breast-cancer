miRNA analysis on breast cancer data
================
Emmanuel Dumont, PhD

## Introduction to the problem

Like many areas in genomics (the study of the genome sequence) and
transcriptomics (the study of the genome expression), we look at a
high-dimensionality problem where the number of features is much larger
than the number of samples in the dataset.

In this particular case, our dataset comes from the expression count of
micro-RNAs (“miRNAs”) that were extracted from formalin-fixed
paraffin-embedded tissue specimens. All tissue specimens come from
women’s breast diagnosed with cancer. “Control” tissues have
non-metastatic breast cancer while “Case” tissues have metastatic breast
cancer. Tissues were paired in (Case, Control) when the cancers had the
same tumor grade and came from women with the age. miRNAs are
single-stranded non-coding RNA with about two dozen nucleotides.

The dataset comprises of 290 samples, 210 of which are paired (samples
came from women with the same age and same tumor type but different
metastatic progression). The goal of this project is to find a set of
miRNAs whose expression levels can be used to differentiate between
healthy patients and patients with breast cancer.

The work is divided in three steps: 1. Prepare the dataset and divide it
into training and test sets. 2. Clustering the miRNAs 3. Apply
machine-learning models on representatives of the clusters.

As we will see below, a visualization method widely used in RNA-seq
(t-distributed stochastic neighbor embedding) was not able to visually
identify clusters between Cases and Controls. Therefore, as expected,
all clustering approaches combined with machine-learning models did not
yield any good predictive results.

## Import, clean the data, and filter the data based on coverage

### Import the matrices of miRNA expression counts per sample

In this file, the NAs were replaced with zeros.

After importing, we concatenate the first two columns (mir\_id and
mir\_seq) to have a unique identifier per miRNA. We then rename each row
with the miRNA new name.

``` r
# minimize all caps (better when using the package DESEq2)
colnames(rawCounts) = tolower(colnames(rawCounts))

# Obtain the sample names by grepping "sRNA
sampleNames <- colnames(rawCounts)[grepl("srna", colnames(rawCounts))]

# Check that every sample is named uniquely
if (length(unique(sampleNames)) != length(sampleNames)){
  print("CHECK. There are duplicate names")
}

# Concatenate the first two columns of the miR dataframe to create 
# a unique ID per isoform
if("mir_seq" %in% colnames(rawCounts)) {
  rawCounts$mir_rna <- paste(rawCounts$mir_id, "_", rawCounts$mir_seq, sep = "")
} else {
  rawCounts$mir_rna = rawCounts$mir_id
}

# Re-construct a new set of columns
rawCounts = rawCounts[, c('mir_rna', sampleNames)]

# Remove the "_count" from the samples' names
names(rawCounts) <- gsub("_count", "", names(rawCounts))

# Use the first column as row names
rownames(rawCounts) <- rawCounts$mir_rna
rawCounts <- rawCounts[,-1]

# Delete variables.
rm(sampleNames)

nbSamplesRaw <- ncol(rawCounts)
nbmiRNAsRaw <- nrow(rawCounts)

head(rawCounts)
```

    ##                                     ol_srna_tmm1_k063g ol_srna_tmm1_k063y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm1_k124g ol_srna_tmm1_k124y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm1_k585g ol_srna_tmm1_k585y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm1_k628g ol_srna_tmm1_k628y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm1_k666g ol_srna_tmm1_k666y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm1_k668g ol_srna_tmm1_k668y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm1_k755g ol_srna_tmm1_k755y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm1_k757g ol_srna_tmm1_k757y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm2_k582g ol_srna_tmm2_k582y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm2_k605g ol_srna_tmm2_k605y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm2_k658g ol_srna_tmm2_k658y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  1                  0
    ##                                     ol_srna_tmm2_k723g ol_srna_tmm2_k723y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm2_k745g ol_srna_tmm2_k745y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  2
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm2_k765g ol_srna_tmm2_k765y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm2_k795g ol_srna_tmm2_k795y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm2_k837g ol_srna_tmm2_k837y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm3_k531g ol_srna_tmm3_k531y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm3_k567g ol_srna_tmm3_k567y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm3_k632g ol_srna_tmm3_k632y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm3_k638g ol_srna_tmm3_k638y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm3_k639g ol_srna_tmm3_k639y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm3_k646g ol_srna_tmm3_k646y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm3_k740g ol_srna_tmm3_k740y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm3_k786g ol_srna_tmm3_k786y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm4_k113g ol_srna_tmm4_k113y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm4_k147g ol_srna_tmm4_k147y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm4_k569g ol_srna_tmm4_k569y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm4_k586g ol_srna_tmm4_k586y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm4_k706g ol_srna_tmm4_k706y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm4_k738g ol_srna_tmm4_k738y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  1
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm4_k774g ol_srna_tmm4_k774y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm4_k808g ol_srna_tmm4_k808y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm5_k090g ol_srna_tmm5_k090y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm5_k534g ol_srna_tmm5_k534y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm5_k572g ol_srna_tmm5_k572y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm5_k587g ol_srna_tmm5_k587y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm5_k631g ol_srna_tmm5_k631y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  1
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm5_k699g ol_srna_tmm5_k699y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm5_k701g ol_srna_tmm5_k701y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm5_k785g ol_srna_tmm5_k785y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm6_k148g ol_srna_tmm6_k148y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm6_k576g ol_srna_tmm6_k576y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm6_k584g ol_srna_tmm6_k584y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm6_k598g ol_srna_tmm6_k598y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm6_k614g ol_srna_tmm6_k614y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm6_k673g ol_srna_tmm6_k673y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm6_k772g ol_srna_tmm6_k772y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm6_k798g ol_srna_tmm6_k798y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm7_k558g ol_srna_tmm7_k558y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm7_k561g ol_srna_tmm7_k561y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm7_k570g ol_srna_tmm7_k570y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm7_k583g ol_srna_tmm7_k583y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  1
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm7_k616g ol_srna_tmm7_k616y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm7_k655g ol_srna_tmm7_k655y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm7_k667g ol_srna_tmm7_k667y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm7_k689g ol_srna_tmm7_k689y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm8_k017g ol_srna_tmm8_k017y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm8_k029g ol_srna_tmm8_k029y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm8_k057g ol_srna_tmm8_k057y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm8_k579g ol_srna_tmm8_k579y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm8_k595g ol_srna_tmm8_k595y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm8_k734g ol_srna_tmm8_k734y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm8_k756g ol_srna_tmm8_k756y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm8_k827g ol_srna_tmm8_k827y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm9_k101g ol_srna_tmm9_k101y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm9_k140g ol_srna_tmm9_k140y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm9_k544g ol_srna_tmm9_k544y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm9_k589g ol_srna_tmm9_k589y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm9_k652g ol_srna_tmm9_k652y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm9_k676g ol_srna_tmm9_k676y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm9_k750g ol_srna_tmm9_k750y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm9_k811g ol_srna_tmm9_k811y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                      0                  0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                     0                  0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                       0                  0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                        0                  0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                   0                  0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                  0                  0
    ##                                     ol_srna_tmm10_k091g ol_srna_tmm10_k091y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm10_k126g ol_srna_tmm10_k126y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm10_k612g ol_srna_tmm10_k612y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm10_k635g ol_srna_tmm10_k635y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm10_k636g ol_srna_tmm10_k636y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm10_k650g ol_srna_tmm10_k650y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm10_k695g ol_srna_tmm10_k695y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm10_k722g ol_srna_tmm10_k722y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm11_k080g ol_srna_tmm11_k080y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm11_k168g ol_srna_tmm11_k168y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm11_k535g ol_srna_tmm11_k535y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm11_k538g ol_srna_tmm11_k538y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm11_k543g ol_srna_tmm11_k543y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm11_k573g ol_srna_tmm11_k573y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm11_k590g ol_srna_tmm11_k590y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm11_k766g ol_srna_tmm11_k766y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm12_k033g ol_srna_tmm12_k033y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm12_k129g ol_srna_tmm12_k129y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm12_k213g ol_srna_tmm12_k213y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm12_k553g ol_srna_tmm12_k553y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm12_k575g ol_srna_tmm12_k575y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm12_k588g ol_srna_tmm12_k588y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm12_k742g ol_srna_tmm12_k742y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm12_k793g ol_srna_tmm12_k793y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm15_k728a ol_srna_tmm15_k728b
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm15_k728c ol_srna_tmm15_k754a
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm15_k754b ol_srna_tmm15_k754c
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm15_k761s ol_srna_tmm15_k761t
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         2                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm15_k761u ol_srna_tmm15_k761v
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm15_k868a ol_srna_tmm15_k868b
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm15_k868c ol_srna_tmm15_k914a
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm15_k914b ol_srna_tmm15_k914c
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm16_k606a ol_srna_tmm16_k606b
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm16_k606c ol_srna_tmm16_k753a
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm16_k753b ol_srna_tmm16_k753c
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm16_k864a ol_srna_tmm16_k864b
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm16_k864c ol_srna_tmm16_k879a
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm16_k879b ol_srna_tmm16_k879b2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                    0
    ##                                     ol_srna_tmm16_k879c ol_srna_tmm16_k967a
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm16_k967b ol_srna_tmm16_k967c
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm17_k604a ol_srna_tmm17_k604b
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm17_k604b2 ol_srna_tmm17_k604c
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                   0
    ##                                     ol_srna_tmm17_k604c2 ol_srna_tmm17_k633a
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                   0
    ##                                     ol_srna_tmm17_k633b ol_srna_tmm17_k633c
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm17_k647g ol_srna_tmm17_k647y
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm17_k880a ol_srna_tmm17_k880b
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm17_k880c ol_srna_tmm17_k979a
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm17_k979b ol_srna_tmm17_k979c
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                       0                   0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                      0                   0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                        0                   0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                         0                   0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                    0                   0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                   0                   0
    ##                                     ol_srna_tmm18_k063g2 ol_srna_tmm18_k063y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm18_k124g2 ol_srna_tmm18_k124y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm18_k585g2 ol_srna_tmm18_k585y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm18_k628g2 ol_srna_tmm18_k628y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm18_k666g2 ol_srna_tmm18_k666y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm18_k668g2 ol_srna_tmm18_k668y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm18_k755g2 ol_srna_tmm18_k755y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm18_k757g2 ol_srna_tmm18_k757y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm19_k091g2 ol_srna_tmm19_k140y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm19_k569y2 ol_srna_tmm19_k583g2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm19_k590y2 ol_srna_tmm19_k606c2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm19_k668g2 ol_srna_tmm19_k673g2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm19_k701g2 ol_srna_tmm19_k721g2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm19_k740y2 ol_srna_tmm19_k782y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm19_k793g2 ol_srna_tmm19_k827y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm19_k837y2 ol_srna_tmm19_k914a2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm20_k582g2 ol_srna_tmm20_k582y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm20_k605g2 ol_srna_tmm20_k605y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm20_k658g2 ol_srna_tmm20_k658y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm20_k723g2 ol_srna_tmm20_k723y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm20_k745g2 ol_srna_tmm20_k745y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm20_k765g2 ol_srna_tmm20_k765y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm20_k795g2 ol_srna_tmm20_k795y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0
    ##                                     ol_srna_tmm20_k837g2 ol_srna_tmm20_k837y2
    ## bhv1-mir-B1_CGGgGTTGGCGGCCGtCGG                        0                    0
    ## bhv1-mir-B1_GCGTTGGCGGgCGaCGGGAA                       0                    0
    ## bhv1-mir-B1_GTCCTCGGCGTgGcCGGC                         0                    0
    ## bhv1-mir-B1_TCCTCGGCGcgGGCGGC                          0                    0
    ## bhv1-mir-B1_as_CAGGCtCTgAATGTCAAAG                     0                    0
    ## bhv1-mir-B1_as_CCcGCCCGTCGCCCGCGCgC                    0                    0

Below we look at the distribution of miRNA types. For instance `hsa`
stands for homo sapiens miRNa. Other classes of miRNAs come from
viruses.

``` r
hsaMirna <- as.numeric(length(which(grepl("hsa", rownames(rawCounts)) == TRUE)))

hsaPerc <- round(100*hsaMirna/length(rownames(rawCounts)), 3)

cat("There are", hsaPerc, "% homosapiens mIRNAs")
```

    ## There are 99.77 % homosapiens mIRNAs

### Import the sample info file

In the matrix above, we do not know if a sample is “Control”
(non-metastatic breast cancer) or “Case” (metastatic breast cancer). The
sample info file will give us this information.

There are 137 “Case” samples and 125 “Control” samples. There is another
sample info file where 224 samples were balanced between case and
control samples and paired with each other. However replicates were
present in the file and they need to be removed. Samples were paired if
they came from women who have the same age and same tumor type but a
different outcome (metastatic vs non-metastatic breast cancer).

``` r
if (datasetName == 'placenta') {
  sampleFileName <- "raw_data/placenta_sample_info.csv"
} else if (datasetName == 'breast') {
  if (pairedSamples == TRUE) {
    sampleFileName <- "raw_data/sample_info_dc_paired_analysis.csv"
  } else {
    sampleFileName <- "raw_data/sample_info_DC.csv"
  }
}

# We import the file created by CD who flagged the samples to remove from the analysis and flagged the samples as "Case" or "Control"
sampleInfo <- read.csv(file = sampleFileName, header = TRUE)

# Remove the samples that neither case or control
sampleInfo = sampleInfo[sampleInfo$status %in%c('Case', 'Control'), ]


if (datasetName == 'breast') {
  # Keep the columns of interest
  sampleInfo = sampleInfo[, c('subsample','status', 'id', 'exclude')]
  
  # Rename the columns
  colnames(sampleInfo) = c('sample','condition', 'id', 'exclude')
  
  # Remove the samples identified to be excluded
  sampleInfo = sampleInfo[is.na(sampleInfo$exclude), ]
  
  # Write id in lower case  
  sampleInfo$id = tolower(sampleInfo$id)


} else if (datasetName == 'placenta') {
  
  # Rename the columns
  colnames(sampleInfo) = c('sample','condition')
  
  # Keep the samples of the method of interest (plasma or PLAP)
  sampleInfo = sampleInfo[grepl(placentaMethod, sampleInfo$sample),]
}
  
# Write names in lower case
sampleInfo$sample = tolower(sampleInfo$sample)

# Rename samples 

if (datasetName == 'breast') {
  #"OL_sRNA_TMM8_k017Y" which is a "Case" becomes "k017_Case"
  sampleInfo$sampleName <- paste(
        str_extract(sampleInfo$id, "k\\d+"), "_",
        sampleInfo$condition, sep = "")
} else if (datasetName == 'placenta') {
  # ol_srna_nick_maternal_plap_h111 becomes h111_plap_Case
  sampleInfo$sampleName <- paste(
        str_extract(sampleInfo$sample, "h\\d+"), "_",
        str_extract(sampleInfo$sample, "p[a-z]+"), "_",
        sampleInfo$condition, sep = "")
}

if (pairedSamples == TRUE) {
  # Remove the replicate samples (at random between 2 replicates)
  #sampleInfo = sampleInfo[!duplicated(sampleInfo$sample),]

  # Remove duplicate samples
  sampleInfo = sampleInfo[!duplicated(sampleInfo$sampleName),]

}

# At this point the dataset should be balanced.
if (dim(sampleInfo[sampleInfo$condition == 'Case', ])[1] == dim(sampleInfo[sampleInfo$condition == 'Control', ])[1]) {
  cat("The dataset is balanced")
} else {
  cat("The dataset is NOT balanced")
}
```

    ## The dataset is balanced

``` r
cat("There are ", length(sampleInfo[sampleInfo$condition == 'Case', ]$sample), "Case samples", "and", length(sampleInfo[sampleInfo$condition == 'Control', ]$sample), "Control samples")
```

    ## There are  105 Case samples and 105 Control samples

``` r
# Keep only the column of the new name and the status
sampleInfo = sampleInfo[, c('sample', 'sampleName', 'condition')]
colnames(sampleInfo) = c('sampleOldName', 'sampleNewName', 'condition')
```

Now, we remove from the matrix of counts all samples that are not in the
sampleInfo file. After that, we’re left with 262 samples in the matrix
of raw counts.

``` r
# Keep the samples of the count matrix that are in the sample info file.
rawCountsId = rawCounts[names(rawCounts) %in% sampleInfo$sampleOldName]

# Identify the samples that are not paired.
idSamples <- colnames(rawCountsId)
allSamples <- colnames(rawCounts)
excludedSamples <- allSamples[!(allSamples %in% idSamples)]
cat("There are", length(excludedSamples), "samples excluded from the dataset")
```

    ## There are 78 samples excluded from the dataset

``` r
# Rename columns of the count matrix.
tmp <- as.data.frame(colnames(rawCountsId))
tmp$newName = apply(tmp, 1, function(x) sampleInfo[sampleInfo$sampleOldName == x, 'sampleNewName'])

colnames(rawCountsId) = tmp$newName

rm(tmp)
```

### Filter the miRNAs by their coverage and prepare a normalized matrix

We remove miRNAs that do not show at least 50 counts in one sample. We
also remove miRNAs where the count is above the 99% percentile after
removing the poorly-expressed miRNAs.

That leaves us with 5k-30k miRNAs (down from 635k), depending on the
parameters.

``` r
# Parameters
minCount <- 50 
maxPerc <- 0.99

samplesAnalysis <- colnames(rawCountsId)

# # We start by making a copy of the dataframe
rawCountsIdCov <- rawCountsId

# Create a column for the max number of counts for a given miRNA
rawCountsIdCov$maxRow = apply(rawCountsIdCov, 1, function(x) max(x))

# Summary of the max number 
summary(rawCountsIdCov$maxRow)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ##       0.0       1.0       1.0      24.1       2.0 1976210.0

``` r
# Filter the matrix of counts out of miRNAs whose max is not at least minCount
rawCountsIdCov = rawCountsIdCov[rawCountsIdCov$maxRow >= minCount, ] 

# Summary of max counts per miRNA after removing the barely-expressed miRNAs
summary(rawCountsIdCov$maxRow)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##      50      73     126    1654     345 1976210

``` r
# Identify the upper bound for the max number of counts
maxCount <- quantile(rawCountsIdCov$maxRow, maxPerc) 

# Filter the matrix of row counts of miRNAs whose max is less than maxCounts
rawCountsIdCov = rawCountsIdCov[rawCountsIdCov$maxRow <= maxCount, ] 

# Summary of max counts per miRNA after removing the over-expressed miRNAs
#summary(rawCountsIdCov$maxRow)
```

We then normalize the matrix of counts per miRNA (between 0 and 1)

``` r
#-------------------------------
# Prepare a normalized matrix (need to convert to numeric first)
countsNorm <- mutate_all(rawCountsIdCov, function(x) as.numeric(as.character(x)))

# We normalize each row (and we need to use the transpose function)
countsNorm = t(apply(
  countsNorm, 1, function(x) round(
                          (x-min(x))/(max(x)-min(x)), 
                          3)))

# Convert the "large matrix" into a dataframe
countsNorm = as.data.frame(countsNorm)

# Remove the column with the maxRow
rawCountsIdCov = rawCountsIdCov[ , samplesAnalysis]
countsNorm = countsNorm[ , samplesAnalysis]
```

## Split the dataset into training (66.66% of data) and test (33.33% of data) sets.

To make sure we avoid over-fitting, we split the dataset into a training
and a test datasets.

``` r
# set seed to ensure reproducible results
set.seed(183)

if (pairedSamples == TRUE) {
  # For paired analysis, we randomize the pairs of samples
  sampledNames <- sample(unique(str_extract(colnames(countsNorm), "k\\d+")))

} else {
  # For unpaired analysis, we randomize all samples
sampledNames <- sample(colnames(countsNorm))
}


# Take 20% of these sample positions at random for the test dataset
sampleIndices <- 1:length(sampledNames)
testIndices <- sample(sampleIndices, trunc(length(sampleIndices)/4))

testSamples <- sampledNames[testIndices]
trainSamples <- sampledNames[-testIndices]

# Create dataframes for testing and training
countsTest = rawCountsIdCov[grepl(paste(testSamples, collapse = "|"), colnames(rawCountsIdCov))]
countsTrain = rawCountsIdCov[grepl(paste(trainSamples, collapse = "|"), colnames(rawCountsIdCov))]

normTest = countsNorm[grepl(paste(testSamples, collapse = "|"), colnames(countsNorm))]
normTrain = countsNorm[grepl(paste(trainSamples, collapse = "|"), colnames(countsNorm))]
```

## Identify representative miRNAs to be used in predictive models

Using boundaries on the number of counts we were able to decrease the
number of miRNAs from 635k to 28k but there are still 100x more features
than samples, which is termed the “curse of dimensionality”.

First we define a function to determine if a column is a Control or Case
based on its name

``` r
# Function to figure out if the sample is a "Control" or a "Case"
findCondition <- function(x) {
  if ( grepl('Case', x) ) {
    answer = "Case"
  } else {
      answer = "Control"
  } 
  return (answer)}
```

### Plot the data in 2 dimensions using t-Distributed Stochastic Neighbor Embedding

t-Distributed Stochastic Neighbor Embedding (tSNE) is a non-linear
dimensionality reduction technique widely used in single cell data
analysis to visualize high-dimensional data in 2 dimensions.

As shown below the miRNAs are unable to cluster the two populations
(Case and Control), leaving very little hope that we would be able to
build a suitable predictive model.

``` r
set.seed(42)

# Create a dataset where rows are sample and columns miRNA
tsneSet <- as.data.frame(t(normTrain))
tsneSet$sample <- rownames(tsneSet)
tsneSet$condition = apply(tsneSet['sample'], 1, findCondition)
tsneSet = tsneSet[, !names(tsneSet) %in% c('sample') ]

# Create a matrix without the column of sample conditions
tsneMatrix <- as.matrix(tsneSet[,1:length(colnames(tsneSet))-1])

# Calculate tsne
tsneOut <- Rtsne(tsneMatrix, pca = FALSE, perplexity = 1, theta = 0) # Run TSNE

# Plot in 2D
plot(tsneOut$Y,col=as.factor(tsneSet$condition), asp=1)
```

![](miRNA_files/figure-gfm/tsne-1.png)<!-- -->

### Method \#1 (naive): T-test on each miRNA between the Cases and Controls

This method is the most intuitive: for each miRNA, we measure if there
is a significant differently-expressed measurement between the Cases and
the Controls.

``` r
ttestDataset = normTrain

# Lists of samples that are "Control"  and  "Case"
controlsTrain = ttestDataset[grepl("Control", colnames(ttestDataset))]
casesTrain = ttestDataset[grepl("Case", colnames(ttestDataset))]

# Create a column with the row names
ttestDataset$mirna = rownames(ttestDataset)

# Compute the  p-value of a t-test for each miRNA
suppressWarnings(ttestDataset$pValue <- apply(ttestDataset, 1,
                              function(x) round(wilcox.test(
                                as.numeric(controlsTrain[x['mirna'], ]),
                                as.numeric(casesTrain[x['mirna'], ]),
                                paired = pairedSamples,
                                alternative = "two.sided")$p.value, 3)))

ttestDataset$controlMean <- apply(ttestDataset, 1,
                              function(x) round(mean( 
                                as.numeric( controlsTrain[x['mirna'], ] ), 
                                na.rm = FALSE ),3))

ttestDataset$caseMean <- apply(ttestDataset, 1,
                              function(x) round(mean( 
                                as.numeric(casesTrain[x['mirna'], ]), 
                                na.rm = FALSE ), 3))

ttestDataset$effectSize <- round(ttestDataset$caseMean / ttestDataset$controlMean, 3)
                              
# Correct for multiple testing using Benjamini & Hochberg 
# criteria (commented because not used)
ttestDataset$adjPValue <- p.adjust(ttestDataset$pValue, method = "BH")

# Delete intermediary files
#rm(controlsTrain, casesTrain)

# Display statistics on the effect size
summary(ttestDataset$effectSize)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   0.000   0.682   0.862     Inf   1.086     Inf       3

``` r
# Display statistics on the p-value
summary(ttestDataset$pValue)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.1730  0.3850  0.4327  0.6810  1.0000

``` r
#----------------------------------------------------
# Remove the miRNA where the effect size is NA
ttestDatasetFilt = ttestDataset[!is.na(ttestDataset$effectSize), ]

#----------------------------------------------------

# Pick min effect size
minEffectSize = 1

# Pick max p-value
maxPValue = 0.05

# Filter based on effect size
ttestDatasetFilt = ttestDatasetFilt[ttestDatasetFilt$effectSize > minEffectSize |
                           ttestDatasetFilt$effectSize < 1/minEffectSize, ]

# Filter based on p-value
ttestDatasetFilt = ttestDatasetFilt[ttestDatasetFilt$pValue < maxPValue, ]

# Create an array of the selected miRNAs
ttestmiRNA <- rownames(ttestDatasetFilt)
cat("There are", length(ttestmiRNA), "miRNA selected by the t-test methodology")
```

    ## There are 532 miRNA selected by the t-test methodology

### Method \#2 (existing bioinformatics package): DESeq2

``` r
# DESeq2 requires the matrix of raw counts as an input
deseqCounts <- countsTrain + 1 # +1 to avoid errors (the matrix cannot contain zeros)

# Sample info file for DESeq2
deseqSampleInfo <- as.data.frame(colnames(deseqCounts))
colnames(deseqSampleInfo) = c('sample')
deseqSampleInfo$condition = apply(deseqSampleInfo['sample'], 1, findCondition)

## DESeq2 Analysis
dds <- DESeqDataSetFromMatrix(deseqCounts, 
                                  colData = deseqSampleInfo, 
                                  design = ~ condition)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
dds$condition <- relevel(dds$condition, ref = "Control")

# Compute the up and down regulated miRNAs
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## -- note: fitType='parametric', but the dispersion trend was not well captured by the
    ##    function: y = a/x + b, and a local regression fit was automatically substituted.
    ##    specify fitType='local' or 'mean' to avoid this message next time.

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 1020 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
## DESeq2 results
deseqModel <- results(dds, filterFun = ihw, alpha = 0.05, name = "condition_Case_vs_Control")
summary(deseqModel)
```

    ## 
    ## out of 8402 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 193, 2.3%
    ## LFC < 0 (down)     : 159, 1.9%
    ## outliers [1]       : 0, 0%
    ## [1] see 'cooksCutoff' argument of ?results
    ## see metadata(res)$ihwResult on hypothesis weighting

``` r
deseqmiRNATable <- as.data.frame(deseqModel)

## Function to grab results
get_upregulated <- function(df){
    key <- intersect(rownames(df)[which(df$log2FoldChange>=1)],
              rownames(df)[which(df$pvalue<=0.05)])

    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    return(results)
  }
get_downregulated <- function(df){
  key <- intersect(rownames(df)[which(df$log2FoldChange<=-1)],
            rownames(df)[which(df$pvalue<=0.05)])

  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

mirnaUpregDESeq <- get_upregulated(deseqmiRNATable)
mirnaDownregDESeq <- get_downregulated(deseqmiRNATable)

# Normalized counts from DESeq matrix
deseqNormTrain <- counts(dds, normalized = T)

deseqmiRNA <- c(rownames(mirnaUpregDESeq), rownames(mirnaDownregDESeq)) 

cat("There are", length(deseqmiRNA), "miRNAs selected by DESeq2")
```

    ## There are 221 miRNAs selected by DESeq2

``` r
#miR_upreg$miRNA_id <- rownames(miR_upreg)
#miR_downreg$miRNA_id <- rownames(miR_downreg)
#miR_upreg <- miR_upreg[,c(8,1,2,3,4,5,6,7)]
#miR_downreg <- miR_downreg[,c(8,1,2,3,4,5,6,7)]

# Write the results in txt files
#dir.create("deseq2results/")
#write.table(deseqNormTrain, "deseq2results/deseqNormTrain.txt", quote = F, sep = "\t")
#write.table(mirnaUpregDESeq, "deseq2results/mirnaUpregDESeq.txt", quote = F, sep = "\t", row.names = F)
#write.table(mirnaDownregDESeq, "deseq2results/mirnaDownregDESeq.txt", quote = F, sep = "\t", row.names = F)
```

### Method \#3: HDBSCAN clustering

HDBSCAN is a clustering technique relying on density (in high
dimensions) widely used in single-cell RNA-seq. When applied to our
dataset, it is unable to find any cluster.

``` r
# Create a dataset where rows are sample and columns miRNA
hdbSet <- as.data.frame(t(normTrain))
hdbSet$sample <- rownames(hdbSet)
hdbSet$condition = apply(hdbSet['sample'], 1, findCondition)
hdbSet = hdbSet[, !names(hdbSet) %in% c('sample') ]

# Create a matrix without the column of sample conditions
hdbSet <- as.matrix(hdbSet[,1:length(colnames(hdbSet))-1])

clusters <- hdbscan(hdbSet, minPts = 10)

# Display the number of clusters -- Could not find any clusters
clusters
```

    ## HDBSCAN clustering for 158 objects.
    ## Parameters: minPts = 10
    ## The clustering contains 0 cluster(s) and 158 noise points.
    ## 
    ##   0 
    ## 158 
    ## 
    ## Available fields: cluster, minPts, cluster_scores, membership_prob,
    ##                   outlier_scores, hc

### Filter datasets with the selected miRNAs

``` r
# Pick the miRNAs to be used in the models: ttestmiRNA, deseqmiRNA
mirnaModel <- deseqmiRNA # deseqmiRNA OR ttestmiRNA OR hdbmiRNA

# Matrices of counts
countsTrainFilt = countsTrain[mirnaModel, ]
countsTestFilt = countsTest[mirnaModel, ]

# Matrices of normalized counts
normTestFilt <- normTest[mirnaModel, ]
normTrainFilt <- normTrain[mirnaModel, ]
```

## Predictions using machine-learning techniques

### Prepare the training and test datasets for the models

``` r
# Transpose the datasets for feeding the models
trainingSet <- as.data.frame(t(normTrainFilt))
testSet <- as.data.frame(t(normTestFilt))

# Create a column with the samples names
trainingSet$sample = rownames(trainingSet)
testSet$sample = rownames(testSet)

# Apply function to create a new column
trainingSet$condition = apply(trainingSet['sample'], 1, findCondition)
testSet$condition = apply(testSet['sample'], 1, findCondition)

# New column names (re-ordered)
trainingSet = trainingSet[, c('condition', mirnaModel)]
testSet = testSet[, c('condition', mirnaModel)]

# Create a vector of the labels in the test dataset (control = 0, case = 1)
testLabels = ifelse(testSet$condition == "Control", 0, 1)
trainingLabels = ifelse(trainingSet$condition == "Control", 0, 1)
```

### Visualize the data using multi-dimensional scaling on the training dataset

``` r
# Create MDS dataset in 2 dimensions
mds <- trainingSet %>%
  dist() %>%
  cmdscale() %>%
  as_tibble()
```

    ## Warning in dist(.): NAs introduced by coercion

    ## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if `.name_repair` is omitted as of tibble 2.0.0.
    ## Using compatibility `.name_repair`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

``` r
colnames(mds) <- c("Dim.1", "Dim.2")
# Add the condition to the dataframe
mds$condition <- trainingSet$condition
mds$condition = as.factor(mds$condition)

# Plot MDS for all data
p <- ggscatter(mds, x = "Dim.1", y = "Dim.2",
        size = 2,
        alpha = 0.5,
        color = 'condition',
        palette =  c("#00AFBB", "#FC4E07"),
        repel = TRUE
        )

print(p)
```

![](miRNA_files/figure-gfm/mds-1.png)<!-- -->

### Lasso logistic regression

Because sequencing of miRNAs is cumbersome and expensive, our goal is to
select as few features as possible. Therefore, we run a penalized
logistic regression using the lasso regression. In this regression, the
coefficients of some less contributive variables are forced to be
exactly zero. Only the most significant variables are kept in the final
model.

``` r
# Dumy code categorical predictor variables
 xTraining <- model.matrix(condition~., trainingSet)[,-1]

# We build cvLasso to find the optimal Lambda
cvLasso <- cv.glmnet(xTraining, trainingLabels, family = "gaussian")

# Build model (Family can be binomial, poisson, gaussian)
lassoModel <- glmnet(xTraining, trainingLabels, alpha = 1, family = "binomial", lambda = cvLasso$lambda.min)

xTest <- model.matrix(condition ~., testSet)[,-1]
probabilities <- lassoModel %>% predict(newx = xTest, type="response")
predictedClasses <- ifelse(probabilities > 0.5, "Case", "Control")

# Confusion matrix
table(pred = predictedClasses, true = testSet[, c('condition')])
```

    ##          true
    ## pred      Case Control
    ##   Case       9      10
    ##   Control   17      16

``` r
# Model accuracy
modelAccuracy = mean(predictedClasses == testSet$condition)

# ROC curve
pred <- prediction(as.vector(probabilities), as.vector(testLabels))
perf <- performance(pred,"tpr","fpr")
par(pty="s")


# Plot the ROC curve
plot(perf,  main = "ROC curve")
```

![](miRNA_files/figure-gfm/lasso-1.png)<!-- -->

``` r
#lines(c(0,1),c(0,1),col = "gray", lty = 4 )

# plot the no-prediction line
lassoAUC <- performance(pred, measure = "auc")
lassoAUC <- lassoAUC@y.values[[1]]

cat("The accuracy is", modelAccuracy, "and the AUC is", lassoAUC)  
```

    ## The accuracy is 0.4807692 and the AUC is 0.5399408

``` r
# Regression parameters  
#coef(cv_lasso, cv_lasso$lambda.min)


# precision/recall curve (x-axis: recall, y-axis: precision)
perf <- performance(pred, "prec", "rec")
#plot(perf)

# sensitivity/specificity curve (x-axis: specificity,
# y-axis: sensitivity)
perf <- performance(pred, "sens", "spec")
#plot(perf)

# Coefficients
#coef(cvLasso, cvLasso$lambda.min)
```

### Support vector machines

``` r
# Treat the condition as a factor
trainingSet$condition = as.factor(trainingSet$condition)

# SVM kernels are polynomial, linear, radial, sigmoid
svmModel <- svm(condition ~ ., data = trainingSet, type = 'C-classification', 
                kernel = "linear", cost = 100, gamma = 1, probability = TRUE)

# Generate predictions on the test dataset
svmPred <- predict(svmModel, 
                    testSet[ , !names(testSet) %in% c('condition')], 
                   probability = TRUE)

# Confusion matrix
table(pred = svmPred, true = testSet[, c('condition')])
```

    ##          true
    ## pred      Case Control
    ##   Case      16       5
    ##   Control   10      21

``` r
# Model accuracy
svmAccuracy = mean(svmPred == testSet$condition)

# ROC
pred <- prediction(as.data.frame(attr(svmPred, "probabilities"))$Case, testLabels)
perf <- performance(pred,"tpr","fpr")
par(pty="s")

# Plot the ROC curve
plot(perf,  main = "ROC curve")
```

![](miRNA_files/figure-gfm/svm-1.png)<!-- -->

``` r
# plot the no-prediction line
#lines(c(0,1),c(0,1),col = "gray", lty = 4 )
svmAUC <- performance(pred, measure = "auc")
svmAUC <- svmAUC@y.values[[1]]

cat("The accuracy is", svmAccuracy,"and the AUC is", svmAUC)  
```

    ## The accuracy is 0.7115385 and the AUC is 0.7115385

``` r
# Regression parameters  
#coef(cv_lasso, cv_lasso$lambda.min)


# precision/recall curve (x-axis: recall, y-axis: precision)
perf <- performance(pred, "prec", "rec")
#plot(perf)

# sensitivity/specificity curve (x-axis: specificity,
# y-axis: sensitivity)
perf <- performance(pred, "sens", "spec")
#plot(perf)
```

### Regression trees

``` r
## Regression tree
rpartModel <- rpart(condition ~ ., data = trainingSet)
rpartPred <- predict(rpartModel, 
                     testSet[ , !names(testSet) %in% c('condition')], type = "class")

# Confusion matrix for rpart
table(pred = rpartPred, true = testSet[, c('condition')])
```

    ##          true
    ## pred      Case Control
    ##   Case      10      12
    ##   Control   16      14

``` r
# Accuracy
rpartAccuracy = mean(rpartPred == testSet$condition)
cat("The accuracy is", rpartAccuracy, "\n")
```

    ## The accuracy is 0.4615385

## Conclusion

The miRNAs were not predictive of the chances of a metastatic breast
cancer.
