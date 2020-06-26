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
gsms <- paste0("11111111111112222332322322332222233232346234623424",
               "42442442242242242252742262688527427742774527452749",
               "33193349334355555555532222222333333336666666117111",
               "55555557888867777777")

gr<- c(rep("AML", 13), rep("Leukemia-B cell", 4), rep("Leukemia-T cell", 2), "Leukemia-B cell", "Leukemia-T cell",
       rep("Leukemia-B cell", 2), "Leukemia-T cell", rep("Leukemia-B cell", 2), rep("Leukemia-T cell", 2), rep("Leukemia-B cell", 5), 
       rep("Leukemia-T cell", 2), "Leukemia-B cell", "Leukemia-T cell", "Leukemia-B cell", "Leukemia-T cell", "AML",
       "Normal-Granulocytes", "Leukemia-B cell", "Leukemia-T cell", "AML", "Normal-Granulocytes", "Leukemia-B cell", "Leukemia-T cell", "AML",
       "Leukemia-B cell", rep("AML", 2), "Leukemia-B cell", rep("AML", 2), "Leukemia-B cell", rep("AML", 2),
       rep("Leukemia-B cell", 2), "AML", rep("Leukemia-B cell", 2), "AML", rep("Leukemia-B cell", 2), "AML", rep("Leukemia-B cell", 2),
       "Normal-B cell", "Leukemia-B cell", "Normal-T cell", "AML", rep("Leukemia-B cell", 2), "Normal-Granulocytes",
       "Leukemia-B cell", "Normal-Granulocytes", rep("Normal-Monocytes", 2), "Normal-B cell", "Leukemia-B cell",
       "Normal-T cell", "AML", "Leukemia-B cell", rep("Normal-T cell", 2), "AML", "Leukemia-B cell", rep("Normal-T cell", 2), 
       "AML", "Normal-B cell", "Leukemia-B cell", "Normal-T cell", "AML", "Normal-B cell", "Leukemia-B cell", "Normal-T cell",
       "AML", "Normal-CD34", rep("Leukemia-T cell", 2), "AML", "Normal-CD34", rep("Leukemia-T cell", 2), "AML", "Normal-CD34",
       rep("Leukemia-T cell", 2), "AML", "Leukemia-T cell", rep("Leukemia-B cell", 9), "Leukemia-T cell", rep("Leukemia-B cell", 7),
       rep("Leukemia-T cell", 8), rep("Normal-Granulocytes", 7), rep("AML", 2), "Normal-T cell", rep("AML", 3), rep("Normal-B cell", 7),
       "Normal-T cell", rep("Normal-Monocytes", 4), "Normal-Granulocytes", rep("Normal-T cell", 7))
       
#extracting expression 
ex <- exprs(gset)

####quality control 

#no log2 scale requierd
#max(ex) min(ex)
#no normalize requierd
pdf("Results/boxplot_subgroups.pdf", width = 64)
boxplot(ex)
dev.off()

#heatmap
pdf("Results/CorHeatmap_subgroups.pdf", width = 20, height = 20)
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
pdf("Results/PCA_samples_subgroups.pdf")
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




