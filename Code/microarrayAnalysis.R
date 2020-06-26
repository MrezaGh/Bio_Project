setwd("D:/Courses/Bio/Project/Project/")

library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(plyr)
library(reshape2)

series <- "GSE48558"

#### load series and platform data from GEO
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
gset <- gset[[1]]

# group names for all samples
gsms <- paste0("1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XXXXX",
               "XXXXXXXXXXXXXXXXXX0X0XXX0X0000X0XX00XX00X0X0X0X0X0",
               "XXX0XXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000110111",
               "00000000000000000000")

gr<- c(rep("AML", 13), rep("Leukemia", 27), "Normal", rep("Leukemia", 3), "Normal", rep("Leukemia", 23),
       "Normal", "Leukemia", "Normal", rep("Leukemia", 3), "Normal", "Leukemia", rep("Normal", 4), "Leukemia", "Normal",
       rep("Leukemia", 2), rep("Normal", 2), rep("Leukemia", 2), rep("Normal", 2), "Leukemia", "Normal", "Leukemia", "Normal",
       "Leukemia", "Normal", "Leukemia", "Normal", "Leukemia", "Normal", rep("Leukemia", 3), "Normal", rep("Leukemia", 3),
       "Normal", rep("Leukemia", 29), rep("Normal", 7), rep("AML", 2), "Normal", rep("AML", 3), rep("Normal", 20))


#extracting expression 
ex <- exprs(gset)

####quality control 

#no log2 scale requierd
#max(ex) min(ex)
#no normalize requierd
pdf("Results/boxplot.pdf", width = 64)
boxplot(ex)
dev.off()

#heatmap
pdf("Results/CorHeatmap.pdf", width = 20, height = 20)
pheatmap(cor(ex), labels_row = gr, labels_col = gr)
dev.off()


####dimention reduction

#principle component Analysis
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

#scaling principal component
ex.scale <- t(scale(t(ex), scale = FALSE))
pc <- prcomp(ex.scale)
pdf("Results/PC_scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()


#principle component for samples
pcr <- data.frame(pc$r[,1:3], Group = gr)
pdf("Results/PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size=3) + theme_bw()
dev.off()


####differential expression analysis

gr <- factor(gr)
gset$description <- gr
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(gr)
fit <- lmFit(gset, design) 
cont.matrix <- makeContrasts(AML-Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("Gene.symbol", "Gene.ID", "adj.P.Val","logFC"))
write.table(tT, "Results/AML_Normal.txt", row.names=F, sep="\t", quote = FALSE)

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file = "Results/AML_Normal_Up.txt", quote = F, row.names = F, col.names = F)
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file = "Results/AML_Normal_Down.txt", quote = F, row.names = F, col.names = F)




