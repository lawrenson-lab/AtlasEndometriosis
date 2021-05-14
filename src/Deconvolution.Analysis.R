
rm(list = ls())

library(MuSiC)
library(xbioc)
library(reshape2)
library(BisqueRNA)
library(biomaRt)
library(UpSetR)
library(Biobase)
library(pheatmap)
library("GEOquery")
require('hgu133plus2')
unloadNamespace("monocle3")
library(cowplot)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'ensembl_transcript_id'), mart = ensembl)

##
title = "25 CCOC"
load("data/ccoc_GSE129617.expr_25_samples.rda")
head(GSE129617.expr)
dim(GSE129617.expr)


df <- data.frame(histotype = c(rep('CCOC', 25)))
rownames(df) <- colnames(GSE129617.expr)

matrix <- GSE129617.expr
matrix <- matrix[!duplicated(rownames(matrix)) & rownames(matrix) %in% 
                   bm$external_gene_name[!is.na(bm$gene_biotype) & 
                                           bm$gene_biotype == 'protein_coding'],]
sd <- apply(X = matrix, MARGIN = 1, FUN = sd)
rank <- (length(sd) - rank(sd) + 1)
select <- rank <= 500

pheatmap(matrix[select,], cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col = df, 
         main = paste0('Top ',sum(select),' protein coding genes (n=', nrow(matrix), ")"))

matrix.eset = ExpressionSet(assayData=matrix, phenoData = new("AnnotatedDataFrame",data=df))

##
GSE9899.anno.clinical <- read.delim("/media/Data01/public_data/CCOC/GSE9899/GSE9899_clinical_anns.csv", sep = ",")
GSE9899.anno.clinical = merge(GSE9899.anno.clinical, read.delim("/media/Data01/public_data/CCOC/GSE9899/AOCSID_GSM_IDs_clinical_anns.csv", sep = ","), by="AOCSID")

title = "GSE9899 - 20 Endometrioid"

GSE9899 <- getGEO('GSE9899',GSEMatrix=TRUE)
GSE9899.expr <- exprs(GSE9899$GSE9899_series_matrix.txt.gz)
head(GSE9899.expr)

annot.BM <-  getBM(mart=ensembl,
                   attributes=c("affy_hg_u133_plus_2", "ensembl_gene_id", "gene_biotype", "external_gene_name"),
                   filter="affy_hg_u133_plus_2",
                   values=rownames(GSE9899.expr),
                   uniqueRows=TRUE)

common = match(rownames(GSE9899.expr), annot.BM$affy_hg_u133_plus_2)
rownames(GSE9899.expr)[!is.na(common)] = annot.BM$external_gene_name[na.omit(common)]

df <- GSE9899.anno.clinical
rownames(df) <- df$GSMID

df <- df[order(df$Primary.Site),]

# Heapmap 0: plot top variable genes
matrix <- GSE9899.expr[, df$GSMID]
matrix <- matrix[complete.cases(matrix),]
matrix <- matrix[!duplicated(rownames(matrix)) & rownames(matrix) %in%
                   bm$external_gene_name[!is.na(bm$gene_biotype) &
                                           bm$gene_biotype == 'protein_coding'],]
sd <- apply(X = matrix, MARGIN = 1, FUN = sd)
rank <- (length(sd) - rank(sd) + 1)
select <- rank <= 50
df = df[,-c(1, 7)]

pheatmap(matrix[select,], cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col = df,
         main = paste0('Top ',sum(select),' protein coding genes (n=', nrow(matrix), ")"))


df = df[df$Subtype == "Endo", , drop=F]
matrix = matrix[, which(colnames(matrix) %in% rownames(df))]

df = df[,3, drop=F]
colnames(df) <- "histotype"

pheatmap(matrix[select,], cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col = df,
         main = paste0('Top ',sum(select),' protein coding genes (n=', nrow(matrix), ")"))


matrix.eset = ExpressionSet(assayData=matrix, phenoData = new("AnnotatedDataFrame",data=df))

######################################

title = "GSE73614 - OvTuCC = 24 OvTuEndo = 35"

five.word <- function(my.string){
  elems = unlist(strsplit(as.character(my.string), " "))
  elems[length(elems)]
}

GSE73614.anno.clinical <- read.delim("/media/Data01/public_data/CCOC/GSE73614/Sample_type.csv", sep = ",")

GSE73614.anno.clinical$histotype <- sapply(GSE73614.anno.clinical$histotype, five.word)
GSE73614.anno.clinical$histotype = gsub("Mayo", "", GSE73614.anno.clinical$histotype)


GSE73614 <- getGEO('GSE73614',GSEMatrix=TRUE)
GSE73614.expr <- exprs(GSE73614$GSE73614_series_matrix.txt.gz)
head(GSE73614.expr)

GSE73614.expr = read.table("/media/Data01/public_data/CCOC/GSE73614/GSE73614_RAW.tar")

annot.BM <-  getBM(mart=ensembl,
                   filters="agilent_wholegenome_4x44k_v2",
                   attributes=c("agilent_wholegenome_4x44k_v2","external_gene_name","ensembl_gene_id","gene_biotype"),
                   values=rownames(GSE73614.expr),
                   uniqueRows=TRUE)
head(annot.BM)

common = match(rownames(GSE73614.expr), annot.BM$agilent_wholegenome_4x44k_v2)
rownames(GSE73614.expr)[!is.na(common)] = annot.BM$external_gene_name[na.omit(common)]
head(GSE73614.expr)


df <- GSE73614.anno.clinical
rownames(df) <- df$SampleName
df = df[,-1, drop=F]

# Heapmap 0: plot top variable genes

matrix <- GSE73614.expr

matrix = matrix[complete.cases(matrix),]
delta1 = max(matrix, na.rm = T) - min(matrix, na.rm = T)
delta2 = max(matrix, na.rm = T) - 0

matrix = (matrix - min(matrix, na.rm = T)) * delta1 / delta2
matrix[1:5, 1:4]

matrix <- matrix[!duplicated(rownames(matrix)) & rownames(matrix) %in% 
                   bm$external_gene_name[!is.na(bm$gene_biotype) & 
                                           bm$gene_biotype == 'protein_coding'],]
sd <- apply(X = matrix, MARGIN = 1, FUN = sd)
rank <- (length(sd) - rank(sd) + 1)
select <- rank <= 50

pheatmap(matrix[select,], cluster_rows=T, show_rownames=F, cluster_cols=F, annotation_col = df, 
         main = paste0('Top ',sum(select),' protein coding genes (n=', nrow(matrix), ")"))


df = df[df$histotype %in% c("OvTuEndo", "OvTuCC"), , drop=F]
matrix = matrix[, which(colnames(matrix) %in% rownames(df))]

pheatmap(matrix[select,], cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col = df, 
         main = paste0('Top ',sum(select),' protein coding genes (n=', nrow(matrix), ")"))


df = df[df$histotype %in% c("OvTuEndo", "OvTuCC"), , drop=F]
df = df[order(df$histotype), , drop=F]
matrix = matrix[, match(rownames(df), colnames(matrix))]

pheatmap(matrix[select,], cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col = df, 
         main = paste0('Top ',sum(select),' protein coding genes (n=', nrow(matrix), ")"))


matrix.eset = ExpressionSet(assayData=matrix, phenoData = new("AnnotatedDataFrame",data=df))



# load
sc = readRDS("rds/epithelial.cells.15.clusters.rds")

levels(sc@meta.data$seurat_clusters) <- c("EnEpi_1", "Meso_1", "EnEpi_2", "EnEpi_3", "UnEpi_1", "Meso_2", "Meso_3", "EnEpi_4",
                                           "EnEpi_5", "Meso_4", "Meso_5", "UnEpi_2", "EnEpi_6", "UnEpi_3", "Meso_6",
                                           "UnEpi_4"
)
Idents(object = sc) <- sc@meta.data$seurat_clusters

prop.cells = data.frame(table(sc@meta.data$active.cluster))
prop.cells <- prop.cells[order(prop.cells$Freq),]
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="Number of cells in each cell type", x ="Cell type", y = "Number of Cells")

prop.cells = data.frame(table(sc@meta.data$Major.Class_2.0))
prop.cells <- prop.cells[order(prop.cells$Freq),]
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="", x ="Major class", y = "Number of Cells")


sc.sel = subset(sc, subset = seurat_clusters %in% c("EnEpi_1", "EnEpi_2", "EnEpi_3", "EnEpi_4", "EnEpi_5", "EnEpi_6"))
sc.sel = subset(sc.sel, subset = Major.Class_2.0 %in% c("Endometrioma", "Eutopic Endometrium", "Extra-ovarian endometriosis"))
single.cell.expression.set <- SeuratToExpressionSet(sc.sel, delimiter='-', position=2, version="v3")
phenoData(single.cell.expression.set) <- AnnotatedDataFrame(sc.sel@meta.data)

prop.cells = data.frame(table(sc.sel@meta.data$active.cluster))
prop.cells <- prop.cells[order(prop.cells$Freq),]
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="Number of cells in each cell type", x ="Cell type", y = "Number of Cells")


prop.cells = data.frame(table(sc.sel@meta.data$Major.Class_2.0))
prop.cells <- prop.cells[order(prop.cells$Freq),]
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="", x ="Major class", y = "Number of Cells")


Est.prop = music_prop(bulk.eset = matrix.eset,
                      sc.eset = single.cell.expression.set, 
                      clusters = 'seurat_clusters',
                      samples = 'SampleName',
                      verbose = T)

bp.esti.pro = Est.prop$Est.prop.weighted
head(bp.esti.pro)

bp.esti.pro = merge(bp.esti.pro, df, by=0)
head(bp.esti.pro)

bp.esti.pro = bp.esti.pro[match(rownames(df), bp.esti.pro$Row.names), ]

collorder = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")


anno_width = unit(3, "cm")
ht_list = rowAnnotation(text = anno_text(bp.esti.pro$Row.names, location = unit(1, "npc"), just = "right", 
                                         gp = gpar(fontsize = 12)))

ht_list = ht_list + rowAnnotation("dist_tss" = anno_barplot(bp.esti.pro[, c(2:7)], bar_width = 1, gp = gpar(fill = collorder, fontsize = 14), 
                                          width = anno_width), show_annotation_name = FALSE) +
  rowAnnotation(Histotype = anno_simple(bp.esti.pro$OvarianHistotype), gp = gpar(fontsize = 20))

draw(ht_list, heatmap_legend_list = Legend(title = "EnEpi Clusters", labels = colnames(bp.esti.pro[, c(2:7)]), legend_gp = gpar(fill = collorder)))






























