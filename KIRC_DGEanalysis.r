#REQUIREMENTS
library(knitr)
library(SummarizedExperiment)
library(edgeR)
library(geneplotter)
library(sva)
#DATA IMPORT
setwd("D:/Documents/BIM/3Q/IEO/PROJECT/IEO")
se <- readRDS(file.path("seKIRC.rds"))
#DGE OBJECT CREATION
dge <- DGEList(counts=assays(se)$counts, genes=mcols(se))
#WITHIN-SAMPLE NORMALIZATION
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
#LOW SEQUENCING DEPTH SAMPLES FILTERING
cpm_cutoff <- 21
mask <- colSums(assays(se)$counts)/1e6 > cpm_cutoff
se.filt <- se[, mask]
dge.filt <- dge[, mask]
#LOWLY EXPRESSED GENES FILTERING
assays(se.filt)$logCPM <- cpm(dge.filt, log=TRUE, prior.count=0.5)
avgexp <- rowMeans(assays(se.filt)$logCPM)
cpm_cutoff <- round(15/min(dge.filt$samples$lib.size/1e+06), digits = 1)
nsamplescutoff <- min(table(se.filt$type))
mask <- rowSums(cpm(dge.filt) > cpm_cutoff) >= nsamplescutoff
se.filt2 <- se.filt[mask, ]
dge.filt2 <- dge.filt[mask, ]
#BETWEEN-SAMPLE NORMALIZATION
dge.filt2 <- calcNormFactors(dge.filt2)
dge.filt2$samples$group <- se.filt2$type
#BATCH EFFECT IDENTIFICATION
tss <- substr(colnames(se.filt2), 6, 7)
center <- substr(colnames(se.filt2), 27, 28)
plate <- substr(colnames(se.filt2), 22, 25)
portionanalyte <- substr(colnames(se.filt2), 18, 20)
samplevial <- substr(colnames(se.filt2), 14, 16)

table(data.frame(TYPE=se.filt2$type, TSS=tss))
table(data.frame(TYPE=se.filt2$type, CENTER=center))
table(data.frame(TYPE=se.filt2$type, PLATE=plate))
table(data.frame(TYPE=se.filt2$type, SAMPLEVIAL=samplevial))
table(data.frame(TYPE=se.filt2$type, PORTIONANALYTE=portionanalyte))

tss_list <- c("A3","B0","B2","B8","CW","CZ")
plate_list <- c(1503,1541,1672)
samplevial_list <- c("01A","11A")
portion <- "01R"
the_mask <- (substr(colnames(dge.filt),6,7) %in% tss_list & substr(colnames(dge.filt),22,25) %in% plate_list & substr(colnames(dge.filt),18,20) == portion & substr(colnames(dge.filt),14,16) %in% samplevial_list)
se.filt3 <- se.filt2[,the_mask]
dge.filt3 <- dge.filt2[,the_mask]
#REMOVAL OF OUTLIERS
se.filt4 <- se.filt3[,grep("5546|5591|4698", colnames(se.filt3), invert=TRUE)]
dge.filt4 <- dge.filt3[,grep("5546|5591|4698", colnames(dge.filt3), invert=TRUE)]
se.dis <- se.filt3[,grep("5546|5591|4698", colnames(se.filt3))]
#SVA
mod <- model.matrix(~ se.filt4$type, colData(se.filt4))
mod0 <- model.matrix(~ 1, colData(se.filt4))
pv <- f.pvalue(assays(se.filt4)$logCPM, mod, mod0)
sv <- sva(assays(se.filt4)$logCPM, mod, mod0)
modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(se.filt4)$logCPM, modsv, mod0sv)
#DEGA
logCPM <- cpm(dge.filt4, log=TRUE, prior.count=3)
tumorExp <- rowMeans(logCPM[, se.filt4$type == "tumor"])
normalExp <- rowMeans(logCPM[, se.filt4$type == "normal"])
log2fc <- tumorExp - normalExp
ranking <- order(abs(log2fc), decreasing = TRUE)

DEG <- data.frame(
  Log2FC = round(log2fc[ranking], digits = 3),
  FC = round(2^log2fc[ranking], digits = 3),
  `1/FC` = round(2^(-log2fc[ranking]), digits = 3),
  logCPM = as.numeric(rowMeans(assays(se.filt4)$logCPM)),
  `p-value` = as.numeric(pvsv), row.names = rowData(se.filt4)$symbol[ranking],
  check.names = FALSE)
DEG$BonfCutoff <- rep(0.05/nrow(DEG), nrow(DEG))
DEG$BonfPvalue <- p.adjust(DEG$`p-value`, method = "bonferroni")
DEG$FDRcutoff <- (1:nrow(DEG) * 0.05)/nrow(DEG)
DEG$FDRpvalue <- p.adjust(DEG$`p-value`, method = "fdr")
plot(DEG$Log2FC, -log10(DEG$`p-value`), pch=".", cex=3, xlab="Log fold-change", ylab="Raw p-value", las=1)
SDEG <- DEG[DEG$FDRpvalue < 1e-4,]
plot(SDEG$Log2FC, -log10(SDEG$FDRpvalue), pch=".", cex=3, xlab="Log fold-change", ylab="Raw p-value", las=1)
SDEG2 <- SDEG[abs(SDEG$Log2FC) > 4 ,]
plot(SDEG2$Log2FC, -log10(SDEG2$FDRpvalue), pch=".", cex=3, xlab="Log fold-change", ylab="Raw p-value", las=1)
abline(h=-log10(max(DEG$`p-value`[SDEG$FDRpvalue <= 0.05])), lty=2)
OEG <- SDEG2[SDEG2$Log2FC > 4,]
UEG <- SDEG2[SDEG2$Log2FC < -4,]