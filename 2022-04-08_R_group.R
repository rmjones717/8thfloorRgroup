#Demo script for 8th floor R group meting 2022-04-08
#based heavily on http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html

#Goal: get from count data to DEG using limma, visualize interactively with Glimma

#First, load in required libs, uncomment lines 8-11 if installs are needed

#if (!requireNamespace("BiocManager"))
#install.packages("BiocManager")
#BiocManager::install(c("limma", "edgeR", "Glimma", "gplots", "RColorBrewer"))
#install.packages(c("tidyverse", "pheatmap"))
#recommend 'usethis' package for working with version control/git


library(tidyverse)
library(edgeR)
library(limma)
library(Glimma)

#read in the data
MabsReads <- read.csv('Mabsc.Expression.Gene.Data.ReadsM.csv')
SampleInfo <- read_tsv('Annotation.txt')

#Rename first column
MabsReads <- MabsReads %>%
  rename("Gene" = X)

#setup as count data and convert to DGEList object - this is a specialized
#list that is used by many tools that determine DGE including limma

countdata <- MabsReads[,-1]
rownames(countdata) <- MabsReads[,1]

y <- DGEList(countdata)

names(y)

y$samples

#add in group from 'Annotation.txt'
group <- factor(SampleInfo$GroupStrain)
y$samples$group <- group

#Convert to CPM and define a threshold
myCPM <- cpm(countdata)
thresh <- myCPM > 0.2
keep <- rowSums(thresh) >=2


plot(myCPM[,1],countdata[,1], ylim=c(0,20), xlim=c(0,3))
abline(v=0.2)

y <- y[keep, keep.lib.sizes = FALSE]

# Get log2 counts per million using the cpm() function from edgeR
logcounts <- cpm(y,log=TRUE)

#normalize using voom method from limma package
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

par(mfrow=c(1,1))
v <- voom(y, design, plot = TRUE)


#look at distributions before and after normalization
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

#generate interactive MDS plot with Glimma

labels <- paste(SampleInfo$SampleID, SampleInfo$GroupStrain)
glMDSPlot(y, labels = labels, group = group, folder="mds")

#set up contrast matrix to specify which comparisons to make for DGE

full_cla_ami_ntp_nca_matrix <- makeContrasts(NTPvCtrl = ntp - ctl, 
                                             AMI_NTPvCtrl = ntp_ami - ctl,
                                             AMIvCtrl = ami - ctl,
                                             CLAvCtrl = cla - ctl,
                                             NTP_CLAvCtrl = ntp_cla - ctl,
                                             NCAvCtrl = nca - ctl,
                                             NTP_AMIvAMI = ntp_ami - ami,
                                             NTP_CLAvCLA = ntp_cla - cla,
                                             levels = design)

fit <- lmFit(v)
fit.full <- contrasts.fit(fit, full_cla_ami_ntp_nca_matrix)
fit.full.eBayes <- eBayes(fit.full)

summary.full <- decideTests(fit.full.eBayes)

summary(summary.full)



#plot interactive volcano plots with Glimma

glXYPlot(x=fit.full.eBayes$coefficients[,1], y=fit.full.eBayes$lods[,1],
         xlab="logFC", ylab="B", main="NTPvCtrl",
         counts=v$E, groups=group, status=summary.full[,1], side.main="Gene", folder="volcano_NTPvCtrl")

glXYPlot(x=fit.full.eBayes$coefficients[,2], y=fit.full.eBayes$lods[,2],
         xlab="logFC", ylab="B", main="NTP_AMIvCtrl",
         counts=v$E, groups=group, status=summary.full[,2],
         side.main="Gene", folder="volcano_NTP_AMIvCtrl")

glXYPlot(x=fit.full.eBayes$coefficients[,3], y=fit.full.eBayes$lods[,3],
         xlab="logFC", ylab="B", main="AMIvCtrl",
         counts=v$E, groups=group, status=summary.full[,3],
         side.main="Gene", folder="volcano_AMIvCtrl")

glXYPlot(x=fit.full.eBayes$coefficients[,4], y=fit.full.eBayes$lods[,4],
         xlab="logFC", ylab="B", main="CLAvCtrl",
         counts=v$E, groups=group, status=summary.full[,4],
         side.main="Gene", folder="volcano_CLAvCtrl")

glXYPlot(x=fit.full.eBayes$coefficients[,5], y=fit.full.eBayes$lods[,5],
         xlab="logFC", ylab="B", main="NTP_CLAvCtrl",
         counts=v$E, groups=group, status=summary.full[,5],
         side.main="Gene", folder="volcano_NTP_CLAvCtrl")

glXYPlot(x=fit.full.eBayes$coefficients[,6], y=fit.full.eBayes$lods[,6],
         xlab="logFC", ylab="B", main="NCAvCtrl",
         counts=v$E, groups=group, status=summary.full[,6],
         side.main="Gene", folder="volcano_NCAvCtrl")

glXYPlot(x=fit.full.eBayes$coefficients[,7], y=fit.full.eBayes$lods[,7],
         xlab="logFC", ylab="B", main="NTP_AMIvAMI",
         counts=v$E, groups=group, status=summary.full[,7],
         side.main="Gene", folder="volcano_NTP_AMIvAMI")

glXYPlot(x=fit.full.eBayes$coefficients[,8], y=fit.full.eBayes$lods[,8],
         xlab="logFC", ylab="B", main="NTP_CLAvCLA",
         counts=v$E, groups=group, status=summary.full[,8],
         side.main="Gene", folder="volcano_NTP_CLAvCLA")
