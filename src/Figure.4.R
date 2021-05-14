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

fibro = readRDS("rds/fibroblasts.cells.20.clusters.rds")
levels(fibro@meta.data$seurat_clusters) <- c("OvF_1", "PerF_1", "PerF_2", "PerF_3", "PerF_4", "EnS_1", "PerF_5", "PerF_6",
                                                                        "PerF_7", "PerF_8", "PerF_9", "EnS_2", "EnS_3", "EnS_4", "PerF_10",
                                                                        "EnS_5", "EnS_6", "PerF_11", "EnF_7", "OvF_2"
)
Idents(object = fibro) <- fibro@meta.data$seurat_clusters

DimPlot(fibro, reduction = "umap", label = T)
DimPlot(fibro, reduction = "umap", label = TRUE, split.by = "Major.Class_2.0", ncol = 5)

##
prop.cells <- data.frame(table(Idents(fibro)))
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="Number of cells in each cluster", x ="Clusters", y = "Number of Cells")

fibro@active.ident <- factor(fibro@active.ident,
                            levels=c("OvF_1", "OvF_2", 
                                     "PerF_1", "PerF_2", "PerF_3", "PerF_4", "PerF_5", "PerF_6", "PerF_7", "PerF_8", "PerF_9", "PerF_10", "PerF_11",
                                     "EnF_1", "EnF_2", "EnF_3", "EnF_4", "EnF_5", "EnF_6", "EnF_7")
)

selected.features = c("MME", "IFITM1", "DCN", "FAP", "COL11A1", "COL1A1", "PDGFRA", "PDGFRB", "ACTA2", "THY1", "CD44", "ESR1", "PGR", "HIF1A")
DotPlot(object = fibro, assay = "RNA", features = selected.features) +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  scale_colour_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b")

##
fibro.markers.df <- FindAllMarkers(fibro, test.use = "MAST")
head(fibro.markers.df, n = 5)

fc.p = 0.5
pv.p = 0.05
fibro.markers.df$STATUS = "NOT.SIG"
fibro.markers.df[fibro.markers.df$avg_logFC < -fc.p & fibro.markers.df$p_val_adj < pv.p, ]$STATUS = "Down"
fibro.markers.df[fibro.markers.df$avg_logFC > fc.p & fibro.markers.df$p_val_adj < pv.p, ]$STATUS = "Up"

##
gs.fibro <- lapply(unique(fibro.markers.df$cluster), function(x) {
  print(x)
  clu <- fibro.markers.df[fibro.markers.df$cluster == x,]
  clu <- clu[clu$p_val_adj < 0.05 & clu$avg_logFC > 0,]
  print(head(clu))
  symbols = clu$gene
  IDs = mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  return(IDs)
})
names(gs.fibro) <- as.character(colnames(data.m.prob.o))
res <- compareCluster(gs.fibro, fun="enrichPathway")
clusterProfiler::dotplot(res, showCategory=3)

# Panel

clus = fibro
index = 1
list.plot = NULL
bsize = 13

##  First column Classes
counts.prop = data.frame(table(clus@meta.data$Major.Class_2.0, clus@meta.data$seurat_clusters))
counts.prop.perc = group_by(counts.prop, Var2) %>% mutate(percent = Freq/sum(Freq))
head(counts.prop.perc)

aux = aggregate(counts.prop$Freq, by=list(Category=counts.prop$Var1), FUN=sum)
aux$Var2 = "Total"
colnames(aux) <- c("Var1", "Freq", "Var2")
aux$percent = aux$Freq / sum(aux$Freq)
counts.prop.perc = rbind(as.data.frame(counts.prop.perc), aux)

list.plot[[index]] <- ggplot(counts.prop.perc, aes(x = Var2, y = percent, fill = Var1)) + 
  geom_bar(position="stack",stat = "identity", width=0.8) +
  scale_fill_manual(name="Class", values = c("Endometrioma" = "#7b3294", 
                                             "Eutopic Endometrium" = "#c2a5cf", 
                                             "Extra-ovarian endometriosis" = "#d9f0d3", 
                                             "No endometriosis detected" = "#a6dba0",
                                             "Unaffected ovary" = "#008837")) +
  scale_x_discrete(limits = rev(levels(counts.prop.perc$Var2))) +
  theme_set(theme_gray(base_size = bsize)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 270, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
index = index + 1

##  Second column Patient
counts.prop = data.frame(table(clus@meta.data$Patient.No._2.0, clus@meta.data$seurat_clusters))
counts.prop.perc = group_by(counts.prop, Var2) %>% mutate(percent = Freq/sum(Freq))
head(counts.prop.perc)

aux = aggregate(counts.prop$Freq, by=list(Category=counts.prop$Var1), FUN=sum)
aux$Var2 = "Total"
colnames(aux) <- c("Var1", "Freq", "Var2")
aux$percent = aux$Freq / sum(aux$Freq)
counts.prop.perc = rbind(as.data.frame(counts.prop.perc), aux)

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
  theme_set(theme_gray(base_size = bsize)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 270, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
index = index + 1

plot_grid(plotlist = list.plot, ncol =  4)

##

fibro
## Session A
class = "peritoneal.stroma"
fibro.sel = subset(fibro, subset = seurat_clusters %in% c("PerF_1", "PerF_2", "PerF_3", "PerF_4", "PerF_5", "PerF_6", "PerF_7", "PerF_8", "PerF_9", "PerF_10", "PerF_11"))
table(fibro.sel@meta.data$seurat_clusters)

## Session B
class = "endometrial.type.stroma"
fibro.sel = subset(fibro, subset = seurat_clusters %in% c("EnS_1", "EnS_2", "EnS_3", "EnS_4", "EnS_5", "EnS_6", "EnS_7"))
table(fibro.sel@meta.data$seurat_clusters)

fibro.sel@meta.data = droplevels(fibro.sel@meta.data)

DimPlot(fibro.sel, reduction = "umap", label = TRUE) +
  ggtitle(class)

# 
Idents(object = fibro.sel) <- fibro.sel@meta.data$Major.Class_2.0

prop.cells <- data.frame(table(Idents(fibro.sel)))
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="Number of cells in each class", x ="Clusters", y = "Number of Cells")

fibro.sel = subset(fibro.sel, subset = Major.Class_2.0 %in% c("Endometrioma", "Eutopic Endometrium", "Extra-ovarian endometriosis"))

#
fibro.sel.DEG <- FindAllMarkers(fibro.sel, test.use = "MAST")
head(fibro.sel.DEG)

fc.p = 0.8
pv.p = 0.05
fibro.sel.DEG$STATUS = "NOT.SIG"
fibro.sel.DEG[fibro.sel.DEG$avg_logFC < -fc.p & fibro.sel.DEG$p_val_adj < pv.p, ]$STATUS = "Down"
fibro.sel.DEG[fibro.sel.DEG$avg_logFC > fc.p & fibro.sel.DEG$p_val_adj < pv.p, ]$STATUS = "Up"

#
gs.fibro.sel <- lapply(unique(fibro.sel.DEG$cluster), function(x) {
  print(x)
  
  clu <- fibro.sel.DEG[fibro.sel.DEG$cluster == x,]
  clu <- clu[clu$p_val_adj < 0.05 & clu$avg_logFC > 0,]
  print(head(clu))
  
  symbols = clu$gene
  IDs = mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  
  return(IDs)
  
})

names(gs.fibro.sel) <- as.character(unique(fibro.sel.DEG$cluster))
res <- compareCluster(gs.fibro.sel, fun="enrichPathway")
clusterProfiler::dotplot(res, showCategory=5, title=class) +
  RotatedAxis()

