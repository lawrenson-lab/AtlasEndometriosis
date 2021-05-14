rm(list = ls())

library(Seurat)
library(ComplexHeatmap)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggrepel)
library('org.Hs.eg.db')
library(ReactomePA)
library(clusterProfiler)
library(circlize)
library("corrplot")

anno <- read.table("files/Endometriosis single cell cohort - Extraovarian.updated.Nov.10.2020.csv", sep = ",", header = TRUE, stringsAsFactors = F)
rownames(anno) <- anno$SampleName
anno

##
epi = readRDS("rds/epithelial.cells.15.clusters.rds")
epi
table(epi@meta.data$Major.Class_2.0)

levels(epi@meta.data$seurat_clusters) <- c("EnEpi_1", "Meso_1", "EnEpi_2", "EnEpi_3", "UnEpi_1", "Meso_2", "Meso_3", "EnEpi_4",
                                           "EnEpi_5", "Meso_4", "Meso_5", "UnEpi_2", "EnEpi_6", "UnEpi_3", "Meso_6",
                                           "UnEpi_4"
)
Idents(object = epi) <- epi@meta.data$seurat_clusters

DimPlot(epi, reduction = "umap", label = TRUE)
DimPlot(epi, reduction = "umap", label = TRUE, split.by = "Major.Class_2.0", ncol = 5)

##

prop.cells <- data.frame(table(Idents(epi)))
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="Number of cells in each cluster", x ="Clusters", y = "Number of Cells")

##
epi@active.ident <- factor(epi@active.ident,
                           levels=c("EnEpi_1", "EnEpi_2", "EnEpi_3", "EnEpi_4", "EnEpi_5", "EnEpi_6", 
                                    "Meso_1", "Meso_2", "Meso_3", "Meso_4", "Meso_5", "Meso_6",
                                    "UnEpi_1", "UnEpi_2", "UnEpi_3", "UnEpi_4")
)
selected.features = c("CDH1", "HIF1A", "LGR5", "EPCAM", "KRT7", "KRT8", "KRT10", "KRT17", "KRT18", "KRT19", "ACTA2","THY1", "CD44", "FOXJ1", "PAX8", "ESR1", "PGR", "WT1", "CALB2", "DES", "PDPN")

DotPlot(object = epi, assay = "RNA", features = selected.features) +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  scale_colour_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b")

##
epi.markers.df <- FindAllMarkers(epi, test.use = "MAST")
head(epi.markers.df, n = 5)

fc.p = 1
pv.p = 0.05

epi.markers.df$STATUS = "NOT.SIG"
epi.markers.df[epi.markers.df$avg_logFC < -fc.p & epi.markers.df$p_val_adj < pv.p, ]$STATUS = "Down"
epi.markers.df[epi.markers.df$avg_logFC > fc.p & epi.markers.df$p_val_adj < pv.p, ]$STATUS = "Up"

gs.ephi <- lapply(unique(epi.markers.df$cluster), function(x) {
  print(x)
  clu <- epi.markers.df[epi.markers.df$cluster == x,]
  clu <- clu[clu$p_val_adj < 0.05 & clu$avg_logFC > 0,]
  print(head(clu))
  symbols = clu$gene
  IDs = mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  return(IDs)
})
names(gs.ephi) <- as.character(colnames(data.m.prob.o))

res <- compareCluster(gs.ephi, fun="enrichPathway")
clusterProfiler::dotplot(res, showCategory=3)

## Panel
clus = epi
index = 1
list.plot = NULL
bsize = 13

# First column Classes
counts.prop = data.frame(table(clus@meta.data$Major.Class_2.0, clus@meta.data$seurat_clusters))
counts.prop.perc = group_by(counts.prop, Var2) %>% mutate(percent = Freq/sum(Freq))
head(counts.prop.perc)

aux = aggregate(counts.prop$Freq, by=list(Category=counts.prop$Var1), FUN=sum)
aux$Var2 = "Total"
colnames(aux) <- c("Var1", "Freq", "Var2")
aux$percent = aux$Freq / sum(aux$Freq)
counts.prop.perc = rbind(as.data.frame(counts.prop.perc), aux)

counts.prop.perc$Var2 = factor(counts.prop.perc$Var2, levels = colnames(data.m.prob.o))

list.plot[[index]] <- ggplot(counts.prop.perc, aes(x = Var2, y = percent, fill = Var1)) + 
  geom_bar(position="stack",stat = "identity", width=0.8) +
  scale_fill_manual(name="Class", values = c("Endometrioma" = "#7b3294", 
                                             "Eutopic Endometrium" = "#c2a5cf", 
                                             "Extra-ovarian endometriosis" = "#d9f0d3", 
                                             "No endometriosis detected" = "#a6dba0",
                                             "Unaffected ovary" = "#008837")) +
  scale_x_discrete(limits = rev(levels(counts.prop.perc$Var2))) +
  coord_flip() +
  theme_set(theme_gray(base_size = bsize)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 270, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
index = index + 1

#  Second column Patient
counts.prop = data.frame(table(clus@meta.data$Patient.No._2.0, clus@meta.data$seurat_clusters))
counts.prop.perc = group_by(counts.prop, Var2) %>% mutate(percent = Freq/sum(Freq))
head(counts.prop.perc)

aux = aggregate(counts.prop$Freq, by=list(Category=counts.prop$Var1), FUN=sum)
aux$Var2 = "Total"
colnames(aux) <- c("Var1", "Freq", "Var2")
aux$percent = aux$Freq / sum(aux$Freq)
counts.prop.perc = rbind(as.data.frame(counts.prop.perc), aux)

counts.prop.perc$Var2 = factor(counts.prop.perc$Var2, levels = colnames(data.m.prob.o))

list.plot[[index]] <- ggplot(counts.prop.perc, aes(x = Var2, y = percent, fill = Var1)) + 
  geom_bar(position="stack",stat = "identity", width=0.8) +
  scale_fill_manual(name="Patient", values = c(  "1" = "#a50026",
                                                 "2" = "#d73027",
                                                 "3" = "#f46d43",
                                                 "4" = "#fdae61",
                                                 "5" = "#fee090",
                                                 "6" = "#ffffbf",
                                                 "7" = "#e0f3f8",
                                                 "8" = "#abd9e9",
                                                 "9"= "#74add1",
                                                 "10" = "#4575b4",
                                                 "11"= "#313695")) +
  scale_x_discrete(limits = rev(levels(counts.prop.perc$Var2))) +
  coord_flip() +
  theme_set(theme_gray(base_size = bsize)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 270, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
index = index + 1

plot_grid(plotlist = list.plot, ncol =  4)


##  Part II  - Endometrial-type Epithelial and Microenvironment analysis

epi
DimPlot(epi, reduction = "umap")
# Session A
class = "endometrial.type.epi"
epi.sel = subset(epi, subset = seurat_clusters %in% c("EnEpi_1", "EnEpi_2", "EnEpi_3", "EnEpi_4", "EnEpi_5", "EnEpi_6"))
epi.sel = subset(epi.sel, subset = Major.Class_2.0 %in% c("Endometrioma", "Eutopic Endometrium", "Extra-ovarian endometriosis"))
table(epi.sel@meta.data$seurat_clusters)

# Session B
class = "microenvironment.epi"
epi.sel = subset(epi, subset = seurat_clusters %in% c("Meso_1", "Meso_2", "Meso_3", "Meso_4", "Meso_5", "Meso_6"))
epi.sel = subset(epi.sel, subset = Major.Class_2.0 %in% c("Endometrioma", "Eutopic Endometrium", "Extra-ovarian endometriosis"))
table(epi.sel@meta.data$seurat_clusters)

DimPlot(epi.sel, reduction = "umap", label = TRUE)

Idents(object = epi.sel) <- epi.sel@meta.data$Major.Class_2.0

##
DimPlot(epi.sel, reduction = "umap", label = TRUE, split.by = "Major.Class_2.0") +
  ggtitle(paste0("Major Class updated ", class))

prop.cells <- data.frame(table(Idents(epi.sel)))
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="Number of cells in each class", x ="Clusters", y = "Number of Cells")

##
epi.sel.DEG <- FindAllMarkers(epi.sel, test.use = "MAST")
head(epi.sel.DEG)

fc.p = 0.8
pv.p = 0.05
epi.sel.DEG$STATUS = "NOT.SIG"
epi.sel.DEG[epi.sel.DEG$avg_logFC < -fc.p & epi.sel.DEG$p_val_adj < pv.p, ]$STATUS = "Down"
epi.sel.DEG[epi.sel.DEG$avg_logFC > fc.p & epi.sel.DEG$p_val_adj < pv.p, ]$STATUS = "Up"

aux.deg = epi.sel.DEG[epi.sel.DEG$STATUS == "Up",]
dot.plot.data = DotPlot(epi.sel, features = unique(aux.deg$gene)) +
  coord_flip() +
  scale_colour_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b") +
  theme(axis.text.x = element_text(angle = 90))
dot.plot.data

features.deg = NULL
for (o in unique(epi.sel.DEG$cluster)) {
  features.deg = c(features.deg, head(epi.sel.DEG$gene[epi.sel.DEG$cluster == o  & epi.sel.DEG$STATUS == "Up"], 10))
}

DotPlot(epi.sel, features = unique(features.deg)) +
  scale_colour_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b") +
  RotatedAxis() +
  coord_flip()

##

gs.epi.sel <- lapply(unique(epi.sel.DEG$cluster), function(x) {
  print(x)
  clu <- epi.sel.DEG[epi.sel.DEG$cluster == x,]
  clu <- clu[clu$p_val_adj < 0.05 & clu$avg_logFC > 0,]
  print(head(clu))
  symbols = clu$gene
  IDs = mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  return(IDs)
})

names(gs.epi.sel) <- as.character(unique(epi.sel.DEG$cluster))

res <- compareCluster(gs.epi.sel, fun="enrichPathway")
n = 7
clusterProfiler::dotplot(res, showCategory=n, title=paste0(class, ": top ", n)) +
  RotatedAxis()

