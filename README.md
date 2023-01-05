# R and Python Investigation: Differential miRNA Expression in the Hippocampus of Mice
## Introduction
Certain probiotic bacteria are known to alter the mood of the animal whose gut they are inhabiting. Lacticaseibacillus rhamnosus JB-1 (JB-1) has previously been shown to alleviate anxiety and depression-like behaviours in stressed mice following feeding, and Limosilactobacillus reuteri 6475 (LR6475) has been shown to help restore social deficits in mice (otherwise referred to as autism-like behaviour). Both of these bacterial activities have been shown to be dependant on Vegus nerve signalling from gut to brain, and in the case of JB-1, immune dependant as well. But with this said, little is known about changes to gene expression in the hippocampus of these mice, the region of the brain associated with emotion, learning, and memory. Previously, in a similar analysis we discovered that there is a great deal of mRNA altered in the blood of these animals prior to feeding with these psychoactive probiotics, but little mRNA is altered in the hippocampus. Here, I analyze nanostring data showing the changes to miRNA expression in the hippocampus of these same mice. For the sake of completeness it should also be noted that blood miRNA was also measured in these mice, but there were no changes of interest between control and experimental groups. 

## Experimental Design
3 groups of 4 mice each were orally gavaged with psychobiotics (JB-1 or LR6475), or phosphate-buffered saline (PBS) control every day for two weeks. 1 hour after the last feeding mice were sacrificed and brains were collected and flash frozen. Hippocapuses were then isolated and lysed to extract total RNA, and total RNA was sent for quality control and miRNA measurement with nanostring using a mouse v1.5 panel. For one of our PBS samples, one of the nanostring spike-in negative controls was abnormally high, and the sample was omitted prior to analysis. Raw .rcc files were loaded into 'R' using the 'NanoStringNorm' package, and background subtraction was performed to remove noise from the count data.

## R Analysis
First, I loaded the raw .rcc files into R
```r
library(NanoStringNorm)
NNS_B2<- newRccSet(
rccFiles = dir("C:/Users/Sandor/Brain Nanostring/RCC files", full.names=TRUE),
rlf = "C:/Users/Sandor/Nanostring Data/RLF/NS_M_miR_v1.5.rlf")
```

Next, I removed background noise from the data
```r
pcn_NNS_B2 <- posCtrlNorm(NNS_B2, summaryFunction = "sum")
bgestimate_B2 <- getBackground(pcn_NNS_B2, bgReference = "negatives", summaryFunction = "mean")
NNS_B2 <- subtractBackground(pcn_NNS_B2, bgEstimates = bgestimate_B2)
```

From this expression set, I then isolated phenotype data and Expression data
```r
B2_expression <- exprs(NNS_B2)
pDataNNS_B2 <- pData(NNS_B2)
```
and exported the CSV files for future analysis in python
```r
write.csv(B2_expression, "C:/Users/Sandor/HippocampalMIRNA_exprs_bg_removed.csv", 
          row.names = T)
write.csv(pDataNNS_B2, "C:/Users/Sandor/HippocampalMIRNA_exprs_bg_removed_pdata.csv", 
          row.names = T)
```
