#Step1 DEGs analysis
library(dplyr)
library(readr)
process_gene_data <- function(file_path, gene_col_name) {
  
  gene_data <- read.delim(file_path, header = TRUE, sep = "\t")
  
  # Filter genes with adjusted P-value < 0.05
  filtered_genes <- gene_data %>% filter(adj.P.Val < 0.05)
  
  # Extract unique gene symbols
  unique_genes <- filtered_genes[[gene_col_name]] %>% unique()
  
  # Remove empty values
  cleaned_genes <- unique_genes[unique_genes != ""]
  
  return(cleaned_genes)
}

# Process each gene data file
ra_genes <- process_gene_data("RA.top.table.tsv", "Gene.Symbol")
ms_genes <- process_gene_data("MS.top.table.tsv", "Gene.symbol")
t1d_genes <- process_gene_data("T1D.top.table.tsv", "Gene.symbol")

# Find the intersection of the gene sets
common_genes <- Reduce(intersect, list(ra_genes, ms_genes, t1d_genes))
print(common_genes)

#Step2 LASSO analysis
library(glmnet)
library(survival)
library(tidyverse)
library(GEOquery)

get_and_preprocess_data <- function(gse_id, gene_column_index, gene_symbol_column) {
  gset <- getGEO(gse_id, destdir=".", AnnotGPL = TRUE, getGPL = TRUE)
  exp <- exprs(gset[[1]])
  gpl <- fData(gset[[1]])
  gpl <- gpl[, c(1, gene_column_index)]
  
  exp <- as.data.frame(exp)
  exp$ID <- rownames(exp)
  exp_symbol <- merge(exp, gpl, by="ID")
  exp_symbol <- na.omit(exp_symbol)
  
  exp_symbol <- exp_symbol[!duplicated(exp_symbol[[gene_symbol_column]]), ]
  datExpr <- t(avereps(exp_symbol[,-c(1, ncol(exp_symbol))], ID=exp_symbol[[gene_symbol_column]]))
  
  return(datExpr)
}

# Define function to perform LASSO regression
perform_lasso <- function(data, common_genes, group, output_prefix) {
  set.seed(1)
  
  # Select common genes
  selected_genes_data <- data[, common_genes, drop = FALSE]
  
  x <- as.matrix(selected_genes_data)
  y <- as.factor(group)
  
  cvla <- glmnet(x, y, family='binomial')
  cv.fit <- cv.glmnet(x, y, family='binomial')
  
  # Save plots
  pdf(paste0(output_prefix, "_LASSO_cvla_plot.pdf"))
  plot(cvla, xvar='lambda', label=TRUE)
  dev.off()
  
  pdf(paste0(output_prefix, "_LASSO_cv_fit_plot.pdf"))
  plot(cv.fit)
  dev.off()
  
  # Extract LASSO gene coefficients
  coef <- coef(cvla, s=cv.fit$lambda.1se)
  index <- which(coef != 0)
  actcoef <- coef[index]
  LassoGene <- rownames(coef)[index]
  
  gene_coef_df <- cbind(Gene=LassoGene, coef=actcoef)
  write.csv(gene_coef_df, paste0(output_prefix, "_LASSO.csv"), row.names=FALSE)
}

# Define main function
main <- function() {
  
  # RA data processing
  ra_data <- get_and_preprocess_data('GSE205962', 15, 'Gene Symbol')
  ra_group <- c(rep(1, 16), rep(0, 4))
  perform_lasso(ra_data, common_genes, ra_group, "RA_GSE205962")
  
  # MS data processing
  ms_data <- get_and_preprocess_data('GSE17048', 3, 'Gene symbol')
  ms_group <- c(rep(0, 45), rep(1, 99))
  perform_lasso(ms_data, common_genes, ms_group, "MS_GSE17048")
  
  # T1D data processing
  t1d_data <- get_and_preprocess_data('GSE44314', 3, 'Gene symbol')
  t1d_group <- c(rep(1, 10), rep(0, 6))
  perform_lasso(t1d_data, common_genes, t1d_group, "T1D_GSE44314")
}

# Run main function
main()


# Step3 RA-single gene GSEA analysis
library(GEOquery)
library(limma)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(enrichplot)
library(dplyr)
RA_gset <- getGEO('GSE205962', destdir=".",AnnotGPL = T,getGPL = T)
RA_exp<-exprs(RA_gset[[1]])
RA_GPL<-fData(RA_gset[[1]])
RA_gpl<- RA_GPL[, c(1, 15)]
RA_exp<-as.data.frame(RA_exp)
RA_exp$ID<-rownames(RA_exp)
RA_exp_symbol<-merge(RA_exp,RA_gpl,by="ID")
RA_exp_symbol<-na.omit(RA_exp_symbol)
table(duplicated(RA_exp_symbol$`Gene Symbol`))
RA_datExpr02<-avereps(RA_exp_symbol[,-c(1,ncol(RA_exp_symbol))],ID=RA_exp_symbol$`Gene Symbol`)
sgene="ROMO1"
data=avereps(RA_datExpr02)
group <- ifelse(data[c(sgene),]>median(data[c(sgene),]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg=topTable(fit2,adjust='fdr',number=nrow(data))
Diff=deg
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file=paste0("1.","RA_DIFF_all.xls"),sep="\t",quote=F,col.names=F)
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(60)){
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=data[afGene,]
Type1=as.data.frame(group)
Type1=Type1[order(Type1$group,decreasing = T),,drop=F]
Type=Type1[,1]
names(Type)=rownames(Type1)
Type=as.data.frame(Type)
anncolor=list(Type=c(High="red",Low="blue"  ))
pheatmap(afExp[,rownames(Type1)],                                                                     
         annotation=Type,                                                           
         color = colorRampPalette(c("blue","white","red"))(50),     
         cluster_cols =F,                                                           
         show_colnames = F,                                                        
         scale="row", 
         fontsize = 10,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=anncolor
)

adjP <- 0.05
aflogFC <- 0.5

Diff$Significant <- ifelse(Diff$P.Value < adjP & Diff$logFC > aflogFC, "Up",
                           ifelse(Diff$P.Value < adjP & Diff$logFC < -aflogFC, "Down", "Not"))

p <- ggplot(Diff, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=Significant), size=4) +  
  scale_color_manual(values=c(pal_npg()(2)[2], "#838B8B", pal_npg()(1))) +
  labs(title = " ") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
  geom_hline(aes(yintercept=-log10(adjP)), colour="gray", linetype="twodash", linewidth=1) + 
  geom_vline(aes(xintercept=aflogFC), colour="gray", linetype="twodash", linewidth=1) +  
  geom_vline(aes(xintercept=-aflogFC), colour="gray", linetype="twodash", linewidth=1)  

p
point.Pvalue=0.05
point.logFc=1.5
logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)
deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))
geneList = data_all_sort$logFC
names(geneList) <- data_all_sort$ENTREZID 
head(geneList)

kk2<- gseKEGG(geneList     = geneList,
              organism = 'hsa',
              nPerm        = 10000,
              minGSSize    = 10,
              maxGSSize    = 200,
              pvalueCutoff = 0.05,
              pAdjustMethod = "none" )
kk1 <- gseGO(geneList = geneList,
             nPerm = 10000,
             minGSSize = 10,
             maxGSSize = 200,
             pvalueCutoff = 0.05,
             pAdjustMethod = "none",
             ont = "BP",
             OrgDb = org.Hs.eg.db)
dotplot(kk1)
class(kk1)
colnames(kk1@result)
kegg_result <- as.data.frame(kk1)
rownames(kk1@result)[head(order(kk1@result$enrichmentScore))]
af1=as.data.frame(kk1@result)
write.table(af1,file=paste0("2.","RA_GSEA_GO.xls"),sep="\t",quote=F,col.names=T)
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file=paste0("2.","RA_GSEA_KEGG.xls"),sep="\t",quote=F,col.names=T)
num=5
gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
#MS & T1D single gene GSEA analysis same to RA

# Step4 RA-Immune Infiltration
library(ggsci)
library(tidyr)
library(ggpubr)
library(stringr) 
library(tidyverse) 
library(pheatmap)
library(ggplot2) 
library(ggpubr)
library(paletteer) 
library(corrplot) 
library(rstatix) 
RA_gset <- getGEO('GSE205962', destdir=".", AnnotGPL = TRUE, getGPL = TRUE)
RA_exp <- exprs(RA_gset[[1]])
RA_GPL <- fData(RA_gset[[1]])
RA_gpl <- RA_GPL[, c(1, 15)]
RA_exp <- as.data.frame(RA_exp)
RA_exp$ID <- rownames(RA_exp)
RA_exp_symbol <- merge(RA_exp, RA_gpl, by="ID")
RA_exp_symbol <- na.omit(RA_exp_symbol)
table(duplicated(RA_exp_symbol$`Gene Symbol`))
RA_datExpr02 <- avereps(RA_exp_symbol[,-c(1, ncol(RA_exp_symbol))], ID=RA_exp_symbol$`Gene Symbol`)
datExpr0 <- t(RA_datExpr02)
b <- read_excel("RA-Immune_Infiltration03.xlsx", sheet = 1)
b <- as.data.frame(b)
rownames(b) <- b[[1]]
b <- b[-1]
a <- read.table("CIBERSORTx_Job1_Results.txt", sep="\t", row.names=1, check.names=FALSE, header=TRUE)
a <- a[,1:22]
identical(rownames(a), rownames(b))
class(b$group)
a$group <- b$group 
a <- a %>% rownames_to_column("sample") 
b <- gather(a, key=CIBERSORT, value=Proportion, -c(group, sample))
my_colors <- c("#66B2FF", "#FFFF00")
ggboxplot(b, x="CIBERSORT", y="Proportion",
          fill="group", palette=my_colors) +
  stat_compare_means(aes(group=group),
                     method="wilcox.test",
                     label="p.signif",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),
                                      symbols=c("***","**","*","ns"))) +
  theme(text=element_text(size=10),
        axis.text.x=element_text(angle=45, hjust=1))
dev.off()
expression_data <- datExpr0[, "ROMO1", drop=FALSE]  
identical(rownames(a), rownames(expression_data))  
resm <- t(a) 
exp_genes <- t(expression_data) 
rb <- rbind(resm, exp_genes) 
rownames(rb) 
rbcor <- cor(t(rb))
rbcor <- rbcor[, !(colnames(rbcor) %in% c("Macrophages M2", "Eosinophils", "sample", "group"))]
rownames(rbcor)
rbcorp <- cor.mtest(rbcor, conf.level=.95)
genes <- "ROMO1" 
p.mat2 <- rbcorp$p
split_cor <- rbcor[genes, , drop=FALSE]
split <- t(split_cor)
splitp <- p.mat2[genes, , drop=FALSE]
splitp <- t(splitp)
mark <- matrix("", nrow = nrow(splitp), ncol = ncol(splitp))
mark[splitp < 0.001] <- "***"
mark[splitp >= 0.001 & splitp < 0.01] <- "**"
mark[splitp >= 0.01 & splitp < 0.05] <- "*"
col2 <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
pheatmap(t(split),
         display_numbers = t(mark),
         number_color = "black",
         fontsize_number = 16,
         color = col2,
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         cluster_rows = FALSE,  
         cluster_cols = FALSE, 
         fontsize_col = 10,
         fontsize_row = 10,
         angle_col = 90)
# MS & T1D Immune Infiltration same to RA
#Validation
# Validation1—Two sample MR
library(TwoSampleMR)
library(pleio)
library(readxl)
library(dplyr)
library(ieugwasr)
ROMO1<- extract_instruments(outcomes = "eqtl-a-ENSG00000106153",p1=5e-06)
dat<-gwasvcf_to_TwoSampleMR(dat)
ROMO1<-dat
filtered_ROMO1 <- subset(dat, pval.exposure< 5e-06)
ROMO1_sorted <- ROMO1[order(ROMO1$P), ]
filtered_ROMO1$`F-statistic` <- (filtered_ROMO1$beta.exposure^2) / (filtered_ROMO1$se.exposure^2)
RA_outcome_dat <- extract_outcome_data(snps = filtered_ROMO1$SNP, outcomes = 'finn-b-RHEUMA_NOS')
MS_outcome_dat <- extract_outcome_data(snps = filtered_ROMO1$SNP, outcomes = 'ieu-a-1024')
T1D_outcome_dat<- extract_outcome_data(snps = filtered_ROMO1$SNP, outcomes = 'ebi-a-GCST90014023')
dat_RA<- harmonise_data(exposure_dat =filtered_ROMO1,  outcome_dat = RA_outcome_dat)
dat_MS<- harmonise_data(exposure_dat =filtered_ROMO1,  outcome_dat = MS_outcome_dat)
dat_T1D<- harmonise_data(exposure_dat =filtered_ROMO1,  outcome_dat = T1D_outcome_dat)
res_RA <- mr(dat_RA)
res_MS <- mr(dat_MS)
res_T1D <- mr(dat_T1D)
pleio_RA<- mr_pleiotropy_test(dat_RA)
OR_RA<-generate_odds_ratios(res_RA)
het_RA <- mr_heterogeneity(dat_RA)
pleio_MS<- mr_pleiotropy_test(dat_MS)
OR_MS<-generate_odds_ratios(res_MS)
het_MS <- mr_heterogeneity(dat_MS)
pleio_T1D<- mr_pleiotropy_test(dat_T1D)
OR_T1D<-generate_odds_ratios(res_T1D)
het_T1D <- mr_heterogeneity(dat_T1D)
# Validation2—External GEO data(Similar to the code provided earlier)
# Validation3—ROSGs data(Similar to the code provided earlier)
# Validation4-Single cell analysis
library(hdf5r)
library(Seurat) 
library(dplyr)
library(multtest)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(future)
library(glmGamPoi)
#RA
data_sample <- Read10X_h5("GSM4819747_RA_filtered_feature_bc_matrix.h5")
data_seurat <- CreateSeuratObject(data_sample, project = "data_sample")
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
#MS
data_seurat<-readRDS("MSpbmc2.rds")
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
#T1D
data_seurat<-readRDS("T1Dpbmc02.rds")
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
#HC
data_seurat<-readRDS("NCpbmc3.rds")
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")

VlnPlot(data_seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
plot1 <- FeatureScatter(data_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
pbmc <- subset(data_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
ncol(GetAssayData(pbmc, assay = "RNA", slot = "counts"))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc))
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dim = 1:20)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:17)
pbmc <- FindClusters(pbmc, resolution = 1.2)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:17)
DimPlot(pbmc, reduction = "umap")
pbmc <- RunTSNE(pbmc, dims = 1:17)
DimPlot(pbmc, reduction = "tsne",label = TRUE)
VlnPlot(pbmc, features = c("ROMO1"), pt.size = 0.001)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
for (i in 0:9) {
  cat("Cluster", i, "markers:\n")
  markers <- FindMarkers(pbmc, ident.1 = i, min.pct = 0.1, logfc.threshold = 0.1)
  print(head(markers))
}
all.genes <- rownames(pbmc)

head(markers)
#RA
# Identify monocytes (cluster 5) - Reference: PMID: 32989127
FeaturePlot(pbmc, features = c("S100A8"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD14"), reduction = "tsne")
FeaturePlot(pbmc, features = c("S100A9"), reduction = "tsne")

# Identify macrophages (cluster 9) - Reference: PMID: 32989127
FeaturePlot(pbmc, features = c("AIF1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("PSAP"), reduction = "tsne")
FeaturePlot(pbmc, features = c("IFITM3"), reduction = "tsne")
FeaturePlot(pbmc, features = c("LST1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("SERPINA1"), reduction = "tsne")

# Identify NK cells (clusters 3, 4, 6) - Reference: PMID: 32989127
FeaturePlot(pbmc, features = c("NCR3"), reduction = "tsne")
FeaturePlot(pbmc, features = c("NKG7"), reduction = "tsne")
FeaturePlot(pbmc, features = c("KLRD1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("GNLY"), reduction = "tsne")

# Identify B cells (clusters 7, 8) - Reference: PMID: 38159454
FeaturePlot(pbmc, features = c("CD79A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD74"), reduction = "tsne")
FeaturePlot(pbmc, features = c("PAX5"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD79B"), reduction = "tsne")

# Identify T cells (clusters 0, 1, 2, 10, 11, 12) - Reference: PMID: 38159454
FeaturePlot(pbmc, features = c("CD8B"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CCR7"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD8A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("IL7R"), reduction = "tsne")
pbmc <- RenameIdents(pbmc, 
                     `0` = "T cells", 
                     `1` = "T cells", 
                     `2` = "T cells",
                     `3` = "NK cells", 
                     `4` = "NK cells", 
                     `5` = "Monocytes",
                     `6` = "NK cells", 
                     `7` = "B cells", 
                     `8` = "B cells",
                     `9` = "Macrophages", 
                     `10` = "T cells", 
                     `11` = "T cells",
                     `12` = "T cells") 
VlnPlot(pbmc, features = c("ROMO1"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("ROMO1"), reduction = "tsne",label = TRUE)
DimPlot(pbmc, reduction = "tsne",label = TRUE)

#MS
# Identify T cells (clusters 0, 1, 3, 4, 5, 8, 9) - Reference: PMID: 38159454
FeaturePlot(pbmc, features = c("CD8B"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CCR7"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD8A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("IL7R"), reduction = "tsne")

# Identify monocytes (cluster 9) - Reference: PMID: 32989127
FeaturePlot(pbmc, features = c("S100A8"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD14"), reduction = "tsne")
FeaturePlot(pbmc, features = c("S100A9"), reduction = "tsne")

# Identify mast cells (cluster 10)
FeaturePlot(pbmc, features = c("GATA2"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CPA3"), reduction = "tsne")

# Identify B cells (cluster 7) - Reference: PMID: 38159454
FeaturePlot(pbmc, features = c("CD79A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD74"), reduction = "tsne")
FeaturePlot(pbmc, features = c("PAX5"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD79B"), reduction = "tsne")

# Identify NK cells (clusters 1, 12, 13, 14) - Reference: PMID: 32989127
FeaturePlot(pbmc, features = c("GZMB"), reduction = "tsne")
FeaturePlot(pbmc, features = c("NKG7"), reduction = "tsne")
FeaturePlot(pbmc, features = c("FCGR3A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("GZMH"), reduction = "tsne")

# Identify macrophages (clusters 0, 15)
FeaturePlot(pbmc, features = c("AIF1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("PSAP"), reduction = "tsne")
FeaturePlot(pbmc, features = c("IFITM3"), reduction = "tsne")
FeaturePlot(pbmc, features = c("LST1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("SERPINA1"), reduction = "tsne")

# Identify erythrocytes (cluster 5)
FeaturePlot(pbmc, features = c("HBA1"), reduction = "tsne")

pbmc <- RenameIdents(pbmc, 
                     `0` = "Macrophages", 
                     `1` = "NK cells", 
                     `2` = "T cells",
                     `3` = "T cells", 
                     `4` = "T cells", 
                     `5` = "RBCs", 
                     `6` = "T cells", 
                     `7` = "B cells", 
                     `8` = "T cells", 
                     `9` = "Monocytes",
                     `10` = "Mast cells", 
                     `11` = "T cells", 
                     `12` = "NK cells", 
                     `13` = "NK cells", 
                     `14` = "NK cells",
                     `15` = "Macrophages")
pbmc <- subset(pbmc, idents = c("RBCs"), invert = TRUE)
VlnPlot(pbmc, features = c("ROMO1"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("ROMO1"), reduction = "tsne",label = TRUE)
DimPlot(pbmc, reduction = "tsne",label = TRUE)

#T1D
# Identify macrophages (cluster 12)
FeaturePlot(pbmc, features = c("AIF1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("PSAP"), reduction = "tsne")
FeaturePlot(pbmc, features = c("IFITM3"), reduction = "tsne")
FeaturePlot(pbmc, features = c("LST1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("SERPINA1"), reduction = "tsne")

# Identify T cells (clusters 0, 1, 2, 4, 6, 11, 15) - Reference: PMID: 38159454
FeaturePlot(pbmc, features = c("CD8B"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CCR7"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD8A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("IL7R"), reduction = "tsne")

# Identify NK cells (cluster 3) - Reference: PMID: 32989127
FeaturePlot(pbmc, features = c("GZMB"), reduction = "tsne")
FeaturePlot(pbmc, features = c("NKG7"), reduction = "tsne")
FeaturePlot(pbmc, features = c("FCGR3A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("GZMH"), reduction = "tsne")

# Identify monocytes (cluster 9) - Reference: PMID: 32989127
FeaturePlot(pbmc, features = c("S100A8"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD14"), reduction = "tsne")
FeaturePlot(pbmc, features = c("S100A9"), reduction = "tsne")

# Identify B cells (cluster 11) - Reference: PMID: 38159454
FeaturePlot(pbmc, features = c("CD79A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD74"), reduction = "tsne")
FeaturePlot(pbmc, features = c("PAX5"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD79B"), reduction = "tsne")

# Identify megakaryocytes (cluster 7) - Reference: PMID: 38159454
FeaturePlot(pbmc, features = c("TUBB1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("GNG11"), reduction = "tsne")
FeaturePlot(pbmc, features = c("GP9"), reduction = "tsne")

# Identify mast cells (cluster 14)
FeaturePlot(pbmc, features = c("GATA2"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CPA3"), reduction = "tsne")

pbmc <- RenameIdents(pbmc, 
                     `0` = "T cells", 
                     `1` = "T cells", 
                     `2` = "T cells",
                     `3` = "NK cells", 
                     `4` = "T cells", 
                     `5` = "NK cells", 
                     `6` = "T cells", 
                     `7` = "NK cells", 
                     `8` = "B cells", 
                     `9` = "Monocytes",
                     `10` = "Megakaryocytes", 
                     `11` = "T cells", 
                     `12` = "Macrophages", 
                     `13` = "NK cells", 
                     `14` = "Mast cells",
                     `15` = "T cells", 
                     `16` = "NK cells")
VlnPlot(pbmc, features = c("ROMO1"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("ROMO1"), reduction = "tsne",label = TRUE)
DimPlot(pbmc, reduction = "tsne",label = TRUE)

#HC
# Identify T cells (clusters 2, 7, 11, 0, 6, 9, 12, 13) - Reference: PMID: 38159454
FeaturePlot(pbmc, features = c("CD8B"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CCR7"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD8A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("IL7R"), reduction = "tsne")

# Identify monocytes (cluster 10) - Reference: PMID: 32989127
FeaturePlot(pbmc, features = c("S100A8"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD14"), reduction = "tsne")
FeaturePlot(pbmc, features = c("S100A9"), reduction = "tsne")

# Identify B cells (cluster 8) - Reference: PMID: 38159454
FeaturePlot(pbmc, features = c("CD79A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD74"), reduction = "tsne")
FeaturePlot(pbmc, features = c("PAX5"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD79B"), reduction = "tsne")

# Identify NK cells (cluster 5) - Reference: PMID: 32989127
FeaturePlot(pbmc, features = c("GZMB"), reduction = "tsne")
FeaturePlot(pbmc, features = c("NKG7"), reduction = "tsne")
FeaturePlot(pbmc, features = c("FCGR3A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("GZMH"), reduction = "tsne")

# Identify macrophages (clusters 1, 3, 4)
FeaturePlot(pbmc, features = c("AIF1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("PSAP"), reduction = "tsne")
FeaturePlot(pbmc, features = c("IFITM3"), reduction = "tsne")
FeaturePlot(pbmc, features = c("LST1"), reduction = "tsne")
FeaturePlot(pbmc, features = c("SERPINA1"), reduction = "tsne")

pbmc <- RenameIdents(pbmc, 
                     `0` = "T cells", 
                     `1` = "Macrophages", 
                     `2` = "T cells",
                     `3` = "Macrophages", 
                     `4` = "Macrophages", 
                     `5` = "NK cells", 
                     `6` = "T cells", 
                     `7` = "T cells", 
                     `8` = "B cells", 
                     `9` = "T cells",
                     `10` = "Monocytes", 
                     `11` = "T cells", 
                     `12` = "T cells", 
                     `13` = "T cells")

VlnPlot(pbmc, features = c("ROMO1"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("ROMO1"), reduction = "tsne",label = TRUE)
DimPlot(pbmc, reduction = "tsne",label = TRUE)
# 3 diseases combined
immune.combined <- readRDS("RA_MS_T1D_pbmc.rds")
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
# Identifying B cells (clusters 6, 7, 16)
FeaturePlot(immune.combined, features = c("CD79B"), reduction = "umap")
FeaturePlot(immune.combined, features = c("CD79A"), reduction = "umap")
FeaturePlot(immune.combined, features = c("CD74"), reduction = "umap")
FeaturePlot(immune.combined, features = c("PAX5"), reduction = "umap")
FeaturePlot(immune.combined, features = c("CD27"), reduction = "umap") # Cluster 7: B cells memory
FeaturePlot(immune.combined, features = c("IL4R"), reduction = "umap") # Clusters 6, 12: B cells naive

# Identifying T cells (clusters 0, 2, 4, 8, 3, 18)
FeaturePlot(immune.combined, features = c("CD8A"), reduction = "umap") # CD8 T cells
FeaturePlot(immune.combined, features = c("CD8B"), reduction = "umap") # CD8 T cells
FeaturePlot(immune.combined, features = c("IL7R"), reduction = "umap")
FeaturePlot(immune.combined, features = c("CCR7"), reduction = "umap") # T cell naive
FeaturePlot(immune.combined, features = c("LTB"), reduction = "umap") # T cell naive
FeaturePlot(immune.combined, features = c("SELL"), reduction = "umap") # T cell naive
FeaturePlot(immune.combined, features = c("CD4"), reduction = "umap") # CD4 T cells
FeaturePlot(immune.combined, features = c("CCR6"), reduction = "umap") # T cell memory
FeaturePlot(immune.combined, features = c("IL17RA"), reduction = "umap") # T cell memory
FeaturePlot(immune.combined, features = c("ICOS"), reduction = "umap")

# Identifying monocytes (clusters 5, 15)
FeaturePlot(immune.combined, features = c("S100A8"), reduction = "umap")
FeaturePlot(immune.combined, features = c("S100A9"), reduction = "umap")
FeaturePlot(immune.combined, features = c("LYZ"), reduction = "umap")

# Identifying NK cells (clusters 1, 3, 14, 20, 19)
FeaturePlot(immune.combined, features = c("GZMB"), reduction = "umap")
FeaturePlot(immune.combined, features = c("GZMH"), reduction = "umap")
FeaturePlot(immune.combined, features = c("NKG7"), reduction = "umap")

# Identifying macrophages (cluster 8: M0, cluster 16: M2)
FeaturePlot(immune.combined, features = c("CD68"), reduction = "umap")
FeaturePlot(immune.combined, features = c("MERTK"), reduction = "umap")
FeaturePlot(immune.combined, features = c("AIF1"), reduction = "umap")
FeaturePlot(immune.combined, features = c("PSAP"), reduction = "umap")
FeaturePlot(immune.combined, features = c("CD80"), reduction = "umap") # M1 macrophages
FeaturePlot(immune.combined, features = c("CD163"), reduction = "umap") # M2 macrophages
FeaturePlot(immune.combined, features = c("NOS2"), reduction = "umap") # M1 macrophages
# Identifying erythrocytes (cluster 11)
FeaturePlot(immune.combined, features = c("HBA1"), reduction = "umap")
# Confirming cluster 9 as Mast cells
FeaturePlot(immune.combined, features = c("CPA3"), reduction = "umap")
FeaturePlot(immune.combined, features = c("GATA2"), reduction = "umap")

# Identifying cluster 20 as Megakaryocytes
FeaturePlot(immune.combined, features = c("GNG11"), reduction = "umap")
FeaturePlot(immune.combined, features = c("TUBB1"), reduction = "umap")

# Confirming clusters 18 and 21 as Plasma cells
FeaturePlot(immune.combined, features = c("IGLL5"), reduction = "umap")
FeaturePlot(immune.combined, features = c("MZB1"), reduction = "umap")
immune.combined1 <- RenameIdents(immune.combined, 
                                 `0` = "T cells CD4 naive", 
                                 `1` = "NK cells", 
                                 `2` = "T cells CD8 naive",
                                 `3` = "NK cells", 
                                 `4` = "T cells CD4 naive", 
                                 `5` = "Monocytes", 
                                 `6` = "B cells naive", 
                                 `7` = "B cells memory", 
                                 `8` = "Macrophages M0", 
                                 `9` = "Mast cells",
                                 `10` = "T cells CD4 naive", 
                                 `11` = "RBCs", 
                                 `12` = "B cells naive", 
                                 `13` = "Mast cells", 
                                 `14` = "NK cells",
                                 `15` = "Monocytes", 
                                 `16` = "Macrophages M2", 
                                 `17` = "T cells CD4 naive",
                                 `18` = "Plasma cells", 
                                 `19` = "Mast cells", 
                                 `20` = "Megakaryocytes", 
                                 `21` = "Plasma cells",
                                 `22` = "Mast cells")

immune.combined1<- subset(immune.combined1, idents = "RBCs", invert = TRUE)
DimPlot(immune.combined1, label = TRUE)
VlnPlot(immune.combined, features = c("ROMO1"))
VlnPlot(immune.combined1, features = c("ROMO1"), pt.size = 0.001)
FeaturePlot(immune.combined, features = c("ROMO1"), reduction = "umap")
#Cell correlation heat map
library(pheatmap)
print(Idents(immune.combined))
immune.combined1@meta.data$celltype <- Idents(immune.combined1)
head(immune.combined1@meta.data$celltype)
Idents(immune.combined1)<- immune.combined1@meta.data$celltype
av.exp<- AverageExpression(immune.combined1)$RNA
features=names(tail(sort(apply(av.exp, 1, sd)),2000))
av.exp<- av.exp[which(row.names(av.exp)%in% features),]
av.exp.mat <- as.matrix(av.exp)
av.exp.cor <- cor(av.exp.mat, method = "spearman")
pheatmap::pheatmap(av.exp.cor)
#Two-sample t-test
# Results from single cell analysis
T1D$group <- 'T1D'
RA$group <- 'RA'
MS$group <- 'MS'
HC$group <- 'HC'
data <- rbind(T1D, RA, MS, HC)
average_HC <- mean(HC$ROMO1, na.rm = TRUE)
average_T1D <- mean(T1D$ROMO1, na.rm = TRUE)
average_RA <- mean(RA$ROMO1, na.rm = TRUE)
average_MS <- mean(MS$ROMO1, na.rm = TRUE)
t_test_T1D <- t.test(T1D$ROMO1, HC$ROMO1, alternative = "two.sided", var.equal = FALSE)
t_test_RA <- t.test(RA$ROMO1, HC$ROMO1, alternative = "two.sided", var.equal = FALSE)
t_test_MS <- t.test(MS$ROMO1, HC$ROMO1, alternative = "two.sided", var.equal = FALSE)




