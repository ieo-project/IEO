#REQUIREMENTS
library(knitr)
library(SummarizedExperiment)
library(edgeR)
library(geneplotter)
library(sva)
library(ggplot2)
#DATA IMPORT
setwd("D:/Documents/BIM/3Q/IEO/PROJECT/IEO/OLD")
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
tss3 <- substr(colnames(se.filt4), 6, 7)
plate3 <- substr(colnames(se.filt4), 22, 25)
portionanalyte3 <- substr(colnames(se.filt4), 18, 20)
samplevial3 <- substr(colnames(se.filt4), 14, 16)
##################DO NOT KNOW WHAT IS THIS USEFUL FOR##########################
#SVA
mod <- model.matrix(~factor(se.filt4$type) + factor(tss3) + factor(plate3), colData(se.filt4)) 
mod0 <- model.matrix(~ 1, colData(se.filt4))
pv <- f.pvalue(assays(se.filt4)$logCPM, mod, mod0)
sv <- sva(assays(se.filt4)$logCPM, mod, mod0)
modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(se.filt4)$logCPM, modsv, mod0sv)
colnames(modsv) <- c(colnames(modsv)[1:9], paste0("SV", 1:sv$n.sv))




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
################################################################################3
mod2 <- model.matrix(~factor(se.filt4$type) + factor(tss3) + factor(plate3), data = colData(se.filt4))
fit <- lmFit(assays(se.filt4)$logCPM,modsv)
fit <- eBayes(fit)
res <- decideTests(fit)
summary(res)
FDRcutoff <- 0.05
res <- decideTests(fit, p.value = FDRcutoff)
summary(res)
tt <- topTable(fit, coef = 2, n = Inf)
head(tt)
genesmd <- data.frame(chr = as.character(seqnames(rowRanges(se.filt4))), symbol = rowData(se.filt4)[,
                                                                                                    1], stringsAsFactors = FALSE, logCPM = as.numeric(rowMeans(assays(se.filt4)$logCPM)))

fit$genes <- genesmd
fit$logCPM <- data.frame(logCPM = as.numeric(rowMeans(assays(se.filt4)$logCPM)), stringsAsFactors = FALSE)

tt <- topTable(fit, coef = 2, n = Inf)
tt.filt <- tt[abs(tt$logFC)>3,]
tt.filt <- tt.filt[tt.filt$P.Value < 1e-5,]
head(tt)

ggplot(data=tt.filt) +
  geom_point(aes(x=logFC,y=-log(P.Value),color=logCPM)) +
  scale_colour_gradientn(colours=c("#000000" ,"#FF0000" )) +
  geom_hline(yintercept=10, linetype="dashed", color = "black",size=1)+
  geom_vline(xintercept = 3, linetype="dashed", 
             color = "black", size=1)+
  geom_vline(xintercept = -3, linetype="dashed", 
             color = "black", size=1)

over <- tt.filt[tt.filt$logFC > 3,]
under <- tt.filt[tt.filt$logFC < -3,]
