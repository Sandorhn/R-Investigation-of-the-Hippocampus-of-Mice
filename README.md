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

Here we more clearly see that if not for the JB-1 outlier all three feed group clusters would be distinct.

It was then time to move on to a differential gene expression analysis, to see if any miRNAs would be differentially expressed between feed groups, or if the insignificant adonis statistic held true.

First, fitting the normalized data to a linear model,
```r
fit_NNS_B2 <- lmFit(NNSB2voom, mmNNS_B2)
```
then defining the feed group contrasts for the linear model, and using an empirical Bayes method to rank genes for differential expression.
```r
contrasts_NNS_B2 <- makeContrasts(JvR = rhamnosus - reuteri,
                                  RvC = reuteri - PBS,
                                  JvC = rhamnosus - PBS,
                                  levels = mmNNS_B2)
contrasts_NNS_B2

fit2_NNS_B2 <- contrasts.fit(fit_NNS_B2, contrasts_NNS_B2)
fit2_NNS_B2 <- eBayes(fit2_NNS_B2)
```
miRNA names were extracted from the original expression set and used to label miRNAs in the linear model.
```r
anno_NNS_B2 <- fData(NNS_B2) 
anno_NNS_B2 <- select(anno_NNS_B2, GeneName)
fit2_NNS_B2$genes <- anno_NNS_B2
```
Then, finally, the differential expression analysis results were examined, assuming a fold change greater than 50% in either direction to be significant.
```r
NNS_results_B2 <- topTable(fit2_NNS_B2, coef = 2, number = Inf)

library(EnhancedVolcano)
Evol_NNS_B2 <- EnhancedVolcano(NNS_results_B2, 
                               lab = NNS_results_B2$GeneNames, 
                               x = 'logFC', 
                               y = 'adj.P.Val',
                               xlim = c(-1.0, 1.1),
                               ylim = c(0, 2.1),
                               xlab = "Log10 Fold of Change", 
                               ylab = "-Log10 adj.P-Value",
                               title = "JB1 vs PBS (adj.P-Value<0.05, FC>1.5X)",
                               pCutoff = 0.05, 
                               FCcutoff = 0.176, 
                               labSize = 4, 
                               labCol = "black",
                               legendLabels = c("Not Sig", 'Sig FC', 'Sig adj.P-Value', "Sig adj.P-Value and FC"),
                               col = c("black", "grey", "grey", "red"), #ooohhh k so it only shows what you define the colours for
                               # drawConnectors = TRUE,
                               # widthConnectors = 0.5
)
Evol_NNS_B2
```
![miRNA Hippo JvC](https://user-images.githubusercontent.com/121974615/210856803-2a2b2afa-d80f-4765-98e8-77a32f6088ed.PNG)
![miRNA Hippo JvR](https://user-images.githubusercontent.com/121974615/210856805-9fa400fe-ae0d-46a0-85f6-2940d1f7d0fa.PNG)
![miRNA Hippo RvC](https://user-images.githubusercontent.com/121974615/210856807-90d6a25f-08b8-4cfb-9b02-2a09990d0b7b.PNG)

But despite distinct PCA feed groups between LR6475 and PBS, no genes were differentially expressed between any groups.

To try to explain the variation between feed groups seen in the PCA, I then conducted a gene set enrichment analysis, to identify if any of the similarly clustered miRNAs were responsible for targeting any of the same mRNA signalling pathways for repression.

In order to do this, I used an R package called RBiomirGS (), which maps miRNAs to their mRNA targets, and performs an enrichment analysis to identify cannonical pathways that are being selectively repressed by the miRNAs.
I conducted seperate analyses for hallmark -
```r
library(RBiomirGS)
mirichB <- select(NNS_results_B2, GeneNames, logFC, adj.P.Val)
rbiomirgs_logistic(objTitle = "NNSGSE_B_hall", mirna_DE = mirichB, 
                   var_mirnaName = "GeneNames", var_mirnaFC = "logFC", var_mirnaP = "adj.P.Val", 
                   mrnalist = Whole, 
                   mrna_Weight = NULL, gs_file = "h.all.v7.4.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = TRUE, clusterType = "PSOCK")
```
and KEGG pathways -
```r
rbiomirgs_logistic(objTitle = "NNSGSE_B_kegg", mirna_DE = mirichB, 
                   var_mirnaName = "GeneNames", var_mirnaFC = "logFC", var_mirnaP = "adj.P.Val", 
                   mrnalist = Whole, 
                   mrna_Weight = NULL, gs_file = "kegg.v5.2.entrez.gmt", optim_method = "IWLS", 
                   p.adj = "fdr", parallelComputing = TRUE, clusterType = "PSOCK")
```
and examined which of those pathways were significantly altered, and in which direction.
```r
rbiomirgs_volcano(gsadfm = NNSGSE_B_kegg_GS, topgsLabel = TRUE, fdr = TRUE, n = 
                    50, gsLabelSize = 3, sigColour = "blue", plotWidth = 500, plotHeight = 
                    500, xLabel = "model coefficient") 
rbiomirgs_volcano(gsadfm = NNSGSE_B_hall_GS, topgsLabel = TRUE, fdr = TRUE, n = 
                    50, gsLabelSize = 3, sigColour = "blue", plotWidth = 500, plotHeight = 
                    500, xLabel = "model coefficient") 
```
![miRNA Hippo JvC Hallmark](https://user-images.githubusercontent.com/121974615/210861053-63285b32-174b-41fa-8f9e-4bb40d854061.PNG)
![miRNA Hippo JvR Hallmark](https://user-images.githubusercontent.com/121974615/210861054-0efa6164-925c-4458-a327-4bfbaf532682.PNG)
![miRNA Hippo RvC Hallmark](https://user-images.githubusercontent.com/121974615/210861055-1c968ab3-96c2-4fa9-bbb3-9c121b67cf46.PNG)
![miRNA Hippo JvC Kegg](https://user-images.githubusercontent.com/121974615/210861086-7c0f303d-012e-4168-8a83-23b0ce746252.PNG)
![miRNA Hippo JvR Kegg](https://user-images.githubusercontent.com/121974615/210861089-531e58e2-8f44-4e84-9302-298653cacec2.PNG)
![miRNA Hippo RvC Kegg](https://user-images.githubusercontent.com/121974615/210861091-fe52c5f0-a379-4d0a-bb99-a98758ce411f.PNG)

Many pathways were found to be significantly altered by miRNA regulation in the hippocampus of mice. Particularly those fed LR6475.
A positive coefficient indicates that there is less miRNA repression of the pathway, while a negative coefficient indicates increased repression of the particular pathways. As you can see, several cannonical signalling pathways involved in inflammation such as the interferon alpha and gamma are altered after psychobiotic feeding, as well as olfactory, following LR6475 feeding, which is interesting because autism is associated with an impaired sense of smell, and LR6475 alleviates autism-like behaviour in mice.

With an explaination of the biological reason for the variance, I then wanted to go back and see if any other machine learning analysis could confirm that these differences between feed groups were real, and for that, I headed over to python for a battery of other machine learning tests.

## Python Analysis



