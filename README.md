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
Another variable, 'Feed' meaning feed group, was then added to the phenotype data
```r
library(stringr)
library(dplyr)
pDataNNS_B2$Feed <- ""
for (i in 1:nrow(pDataNNS_B2)){
  if(str_detect(pDataNNS_B2$SampleID[i], "C"))
  {pDataNNS_B2$Feed[i] <- "PBS"}
  
  if(str_detect(pDataNNS_B2$SampleID[i], "J"))
  {pDataNNS_B2$Feed[i] <- "rhamnosus"}
  
  if(str_detect(pDataNNS_B2$SampleID[i], "R"))
  {pDataNNS_B2$Feed[i] <- "reuteri"}
}
pDataNNS_B2$Feed
```
I then created the model matrix
```r
library(edgeR)
library(limma)
mmNNS_B2 <- model.matrix(~0 + pDataNNS_B2$Feed)
colnames(mmNNS_B2) <- c("PBS", "reuteri", "rhamnosus")
```
and filtered out low-expressed miRNAs from the data.
```r
dgeTMM_B2 <- DGEList(B2_expression, group = pDataNNS_B2$Feed)
keep.exprsB2 <- filterByExpr(dgeTMM_B2, group=pDataNNS_B2$Feed)#ohhh so dont use the logCPM for this
fdgeTMM_B2 <- dgeTMM_B2[keep.exprsB2,, keep.lib.sizes=FALSE]
```
Next, I wanted to normalize the data prior to conducting a PCA.
Nanostring data is treated like count data, and so the appropriate normalization method is the trimmed mean of M-values (TMM) normalization.
```r
fdgeTMM_B2 <- calcNormFactors(fdgeTMM_B2, method = "TMM")
NNSB2voom <- voom(fdgeTMM_B2, mmNNS_B2, plot = TRUE)
```
![Voom Plot](https://user-images.githubusercontent.com/121974615/210846615-8ff07bdc-9730-4d13-8b32-4928466d64e5.png)

The mean-variance trend plot suggests that low-expressed genes were adequately removed, as lower count size is not positively correlated with a higher standard deviation.
Next, I performed a principal componant analysis (PCA)
```r
library(dplyr)
library(ggplot2)
pca_NNS_B2 <- prcomp(t(NNSB2voom$E))
cbind(pDataNNS_B2, pca_NNS_B2$x) %>%
  ggplot(aes(x = PC1, y = PC2, col=Feed, label=paste("", Feed))) + geom_point(aes(size = 4)) + theme_bw(base_size = 15) +
  xlab("PC1: 31% variance") + ylab("PC2: 18% variance") + guides(colour = guide_legend(override.aes = list(size=5)))
```
![PCA](https://user-images.githubusercontent.com/121974615/210848496-3199f948-40f5-46cf-b6a9-ba90b05f4ebb.png)

The PCA suggests that there may be differences in miRNA expression in the hippocampus of PBS and LR6475 (reuteri), but an outlier from the rhamnosus (JB-1) group brings into question if miRNAs in the hippocampuses of the mice it was fed to will be significantly different from the other feed groups.
To gain a greater understanding of the significance of the PCA plot, I conducted a permutational multivariate analysis of variance using adonis.
```r
library(vegan)
hippo.pca.miRNA.data <- t(NNSB2voom$E)
adonis2(hippo.pca.miRNA.data ~ Feed, data = hippo_pca, method = 'eu', perm=100000)
```
This produced an F statistic of 0.081, suggesting that the PCA clusters have a 92% chance of being truely distinct. ALTHOUGH THIS IS CONSIDERED STATISTICALLY INSIGNIFICANT, I CONTINUED THE ANALYSIS FOR THE EXERCISE, and because of personal curiousity, knowing that everywhere other than in science, a 92% chance is pretty good odds. It should be noted as a disclaimer however that it was because of this statistic that this data could not be published. Professionally, I agree with the decision not to publish it, and it was why I decided to use it for this fun R and python exercise here.

I also conducted a 3D PCA to more clearly illustrate the near perfect distinction between each of the three feed groups. 
```r
library(pca3d)
pca3d(pca_NNS_B2, group = pDataNNS_B2$Feed, show.ellipses = F, radius = 5, legend = "topleft")
```
![3d PCA](https://user-images.githubusercontent.com/121974615/210852532-7673c8bd-8295-46c6-a4ef-75e20a3a70bf.PNG)

Here we more clearly see that if not for the JB-1 outlier all three feed group clusters would ben distinct.

It was then time to move on to a differential gene expression analysis, to see if any miRNAs would be differentially expressed between feed groups, or if the insignificant adonis statistic held true.



