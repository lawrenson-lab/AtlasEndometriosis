rm(list = ls())

library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)

anno <- read.table("files/Endometriosis single cell cohort - Extraovarian.updated.Nov.10.2020.csv", sep = ",", header = TRUE, stringsAsFactors = F)
rownames(anno) <- anno$SampleName
anno

############

# Epithelial
sc = readRDS("rds/epithelial.cells.15.clusters.rds")
levels(sc@meta.data$seurat_clusters) <- c("EnEpi_1", "Meso_1", "EnEpi_2", "EnEpi_3", "UnEpi_1", "Meso_2", "Meso_3", "EnEpi_4",
                                          "EnEpi_5", "Meso_4", "Meso_5", "UnEpi_2", "EnEpi_6", "UnEpi_3", "Meso_6",
                                          "UnEpi_4"
)
Idents(object = sc) <- sc@meta.data$seurat_clusters

DimPlot(sc, reduction = "umap", label = TRUE)

## Session A
class = "endometrial.type.epi"
sc.sel = subset(sc, subset = seurat_clusters %in% c("EnEpi_1", "EnEpi_2", "EnEpi_3", "EnEpi_4", "EnEpi_5", "EnEpi_6"))
sc.sel = subset(sc.sel, subset = Major.Class_2.0 %in% c("Endometrioma", "Eutopic Endometrium", "Extra-ovarian endometriosis", "No endometriosis detected"))
prop.cells <- data.frame(table(sc.sel@meta.data$Major.Class_2.0))
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="Number of cells in each Class", x ="Class", y = "Number of Cells")


## Session B
class = "microenvironment.epi.rm11"
sc.sel = subset(sc, subset = seurat_clusters %in% c("Meso_1", "Meso_2", "Meso_3", "Meso_4", "Meso_5", "Meso_6"))
sc.sel = subset(sc.sel, subset = Major.Class_2.0 %in% c("Endometrioma", "Extra-ovarian endometriosis", "Unaffected ovary", "No endometriosis detected"))
prop.cells <- data.frame(table(sc.sel@meta.data$Major.Class_2.0))
prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)

ggplot(prop.cells, aes(Var1, Freq)) +
  geom_col() +
  theme_minimal(base_size = 15) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(title="Number of cells in each Class", x ="Class", y = "Number of Cells")



#################
# Fibroblast
sc = readRDS("rds/fibroblasts.cells.20.clusters.rds")

levels(sc@meta.data$seurat_clusters) <- c("OvF_1", "PerF_1", "PerF_2", "PerF_3", "PerF_4", "EnS_1", "PerF_5", "PerF_6",
                                             "PerF_7", "PerF_8", "PerF_9", "EnS_2", "EnS_3", "EnS_4", "PerF_10",
                                             "EnS_5", "EnS_6", "PerF_11", "EnS_7", "OvF_2"
)
Idents(object = sc) <- sc@meta.data$seurat_clusters

class = "endometrial.type.stroma"
sc.sel = subset(sc, subset = seurat_clusters %in% c("EnS_1", "EnS_2", "EnS_3", "EnS_4", "EnS_5", "EnS_6", "EnS_7"))
sc.sel = subset(sc.sel, subset = Major.Class_2.0 %in% c("Endometrioma", "Eutopic Endometrium", "Extra-ovarian endometriosis", "No endometriosis detected"))

DimPlot(sc.sel, reduction = "umap", label = TRUE)

expression_matrix <- sc.sel@assays[["RNA"]]@counts
dim(expression_matrix)
pd <- sc.sel@meta.data
fData <- data.frame(gene_short_name = row.names(expression_matrix), row.names = row.names(expression_matrix))

# Construct monocle cds
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = pd,
                         gene_metadata = fData)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)

plot_cells(cds,
           label_groups_by_cluster=F,
           label_branch_points=FALSE,
           label_leaves=FALSE,
           color_cells_by = "seurat_clusters",
           group_label_size = 0,
           cell_size = .7,
           rasterize = 5) +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds)

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds, reduction_method = "UMAP")

plot_cells(cds,
           label_groups_by_cluster=F,
           color_cells_by = "Major.Class",
           group_label_size = 0,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = .7,
           rasterize = 5) +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(name="Class", values = c("Endometrioma" = "#7b3294", 
                                             "Eutopic Endometrium" = "#c2a5cf", 
                                             "Extra-ovarian endometriosis" = "#d9f0d3", 
                                             "No endometriosis detected" = "#a6dba0",
                                             "Unaffected ovary" = "#008837")) 

## Step 4: Cluster the cells
cds <- cluster_cells(cds)
plot_cells(cds,
           color_cells_by = "partition",
           group_label_size = 0,
           cell_size = .7,
           rasterize = 5) +
  theme_bw(base_size = 20)

## Step 5: Learn a graph
cds <- learn_graph(cds, use_partition = F)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=T,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size = 7,
           cell_size = .7,
           rasterize = 5) +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

# a helper function to identify the root principal points:
get_earliest_principal_node_cluster <- function(cds, time_bin=c("1", "2")){
  cell_ids <- which(cds@clusters[["UMAP"]]$clusters %in% time_bin)
  print(length(cell_ids))

  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

### host markers
selected.features = c("CD44", "LGR5", "THY1", "PGR", "IFTM1")

DotPlot(object = sc.sel, assay = "RNA", features = selected.features) +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  scale_colour_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b")


# Epi endometrial-type
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node_cluster(cds, time_bin = "3"))

# Epi mesothelial
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node_cluster(cds, time_bin = "3"))

# Fibro endometrial-type
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node_cluster(cds, time_bin = "4"))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = .7,
           graph_label_size=1.5) +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())


# Endometrial-type Epithelial markers
markers = c("MSLN", "KRT8", "SLPI", "LUM")

# Endometrial-type Stroma markers
markers = c("IFITM1", "ESR1", "PGR", "C3", "C7", "STAR", "PLA2G2A", "CXCL2", "CXCL12","ECM1", "CD44")

plot_cells(cds, genes=markers,
           show_trajectory_graph=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = .5,
           rasterize = 2, 
           scale_to_range = T) +
  theme_bw(base_size = 20) +
  scale_colour_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", na.value = "white") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_rect(color="black", fill="white", size=0, linetype="solid"))


