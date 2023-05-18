library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(cowplot)


commit: commit 3384307eda0bb54e237b9cd142490bd5da678778

#create seurat obj
for (r1 in c("wt1","wt3","utxtcd1","utxtcd2")){
  seurat_data1 <- Read10X(data.dir = paste0("../scRNAseq/lisa/round1/", r1))
  seurat_obj1 <- CreateSeuratObject(counts = seurat_data1, 
                                    min.cells = 20,
                                    min.features=200,
                                    project = r1)
  
  assign(r1, seurat_obj1)
}

head(wt1@meta.data)

for (r2 in c("utx596","utx600","wt595","wt599")){
  seurat_data2 <- Read10X(data.dir = paste0("../scRNAseq/lisa/round2/", r2))
  seurat_obj2 <- CreateSeuratObject(counts = seurat_data2, 
                                    min.cells = 20,
                                    min.features=200,
                                    project = r2)
  
  assign(r2, seurat_obj2)
}

head(utx596@meta.data)

for (r3 in c("utx8wk","utx13wk","wt8wk","wt13wk")){
  seurat_data3 <- Read10X(data.dir = paste0("../scRNAseq/lisa/round3/", r3))
  seurat_obj3 <- CreateSeuratObject(counts = seurat_data3, 
                                    min.cells = 20,
                                    min.features=200,
                                    project = r3)
  
  assign(r3, seurat_obj3)
}

head(utx8wk@meta.data)


# Create a merged Seurat object
wt8 <- merge(x = wt1,
             y = c(wt3,wt8wk),
             add.cell.id = c("wt1","wt3","wt8wk")
)

wt13 <- merge(x = wt595,
              y = c(wt599,wt13wk),
              add.cell.id = c("wt595","wt599","wt13wk")
)

utx8 <- merge(x = utxtcd1,
              y = c(utxtcd2,utx8wk), 
              add.cell.id = c("utxtcd1","utxtcd2","utx8wk")
)

utx13 <- merge(x = utx596,
               y = c(utx600,utx13wk), 
               add.cell.id = c("utx596","utx600","utx13wk")
)

# Add number of genes per UMI for each cell to metadata ###########
wt8$log10GenesPerUMI <- log10(wt8$nFeature_RNA) / log10(wt8$nCount_RNA)
wt13$log10GenesPerUMI <- log10(wt13$nFeature_RNA) / log10(wt13$nCount_RNA)
utx8$log10GenesPerUMI <- log10(utx8$nFeature_RNA) / log10(utx8$nCount_RNA)
utx13$log10GenesPerUMI <- log10(utx13$nFeature_RNA) / log10(utx13$nCount_RNA)

# Compute percent mito ratio
wt8$mitoRatio <- PercentageFeatureSet(object = wt8, pattern = "^mt-")
wt8$mitoRatio <- wt8@meta.data$mitoRatio / 100

wt13$mitoRatio <- PercentageFeatureSet(object = wt13, pattern = "^mt-")
wt13$mitoRatio <- wt13@meta.data$mitoRatio / 100

utx8$mitoRatio <- PercentageFeatureSet(object = utx8, pattern = "^mt-")
utx8$mitoRatio <- utx8@meta.data$mitoRatio / 100

utx13$mitoRatio <- PercentageFeatureSet(object = utx13, pattern = "^mt-")
utx13$mitoRatio <- utx13@meta.data$mitoRatio / 100

##################################    WT8 #####
# Create metadata dataframe for wt8
metadata.wt8 <- wt8@meta.data

# Add cell IDs to metadata
metadata.wt8$cells <- rownames(metadata.wt8)

# Rename columns
metadata.wt8 <- metadata.wt8 %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata.wt8$sample <- NA
metadata.wt8$sample[which(str_detect(metadata.wt8$cells, "^wt1_"))] <- "wt1"
metadata.wt8$sample[which(str_detect(metadata.wt8$cells, "^wt3_"))] <- "wt3"
metadata.wt8$sample[which(str_detect(metadata.wt8$cells, "^wt8wk_"))] <- "wt8wk"


metadata.wt8$mice <- NA
metadata.wt8$mice <- "WT"
metadata.wt8$week <- NA
metadata.wt8$week <- "Wk8_WT"

# Add metadata back to Seurat object
wt8@meta.data <- metadata.wt8
head(wt8@meta.data)


metadata.wt8 %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata.wt8 %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata.wt8 %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata.wt8 %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata.wt8 %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata.wt8 %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata.wt8 %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Filter out low quality reads using selected thresholds - these will change with experiment
wt8 <- subset(x = wt8, 
              subset= (nUMI >= 1000) & 
                (nGene >= 500) & 
                (log10GenesPerUMI > 0.8) & 
                (mitoRatio < 0.20))

##################################   WT13 ####
# Create metadata dataframe for wt8
metadata.wt13 <- wt13@meta.data

# Add cell IDs to metadata
metadata.wt13$cells <- rownames(metadata.wt13)

# Rename columns
metadata.wt13 <- metadata.wt13 %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata.wt13$sample <- NA
metadata.wt13$sample[which(str_detect(metadata.wt13$cells, "^wt595_"))] <- "wt595"
metadata.wt13$sample[which(str_detect(metadata.wt13$cells, "^wt599_"))] <- "wt599"
metadata.wt13$sample[which(str_detect(metadata.wt13$cells, "^wt13wk_"))] <- "wt13wk"

metadata.wt13$mice <- NA
metadata.wt13$mice <- "WT"

metadata.wt13$week <- NA
metadata.wt13$week <- "Wk13_WT"

# Add metadata back to Seurat object
wt13@meta.data <- metadata.wt13

metadata.wt13 %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata.wt13 %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata.wt13 %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata.wt13 %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata.wt13 %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata.wt13 %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata.wt13 %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Filter out low quality reads using selected thresholds - these will change with experiment
wt13 <- subset(x = wt13, 
               subset= (nUMI >= 1100) & 
                 (nGene >= 600) & 
                 (log10GenesPerUMI > 0.8) & 
                 (mitoRatio < 0.20))

##################################   UTX8 #########
# Create metadata dataframe for wt8
metadata.utx8 <- utx8@meta.data

# Add cell IDs to metadata
metadata.utx8$cells <- rownames(metadata.utx8)

# Rename columns
metadata.utx8 <- metadata.utx8 %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata.utx8$sample <- NA
metadata.utx8$sample[which(str_detect(metadata.utx8$cells, "^utxtcd1_"))] <- "utxtcd1"
metadata.utx8$sample[which(str_detect(metadata.utx8$cells, "^utxtcd2_"))] <- "utxtcd2"
metadata.utx8$sample[which(str_detect(metadata.utx8$cells, "^utx8wk_"))] <- "utx8wk"

metadata.utx8$mice <- NA
metadata.utx8$mice <- "UTX"

metadata.utx8$week <- NA
metadata.utx8$week <- "Wk8_UTX"

# Add metadata back to Seurat object
utx8@meta.data <- metadata.utx8

metadata.utx8 %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata.utx8 %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata.utx8 %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata.utx8 %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata.utx8 %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata.utx8 %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata.utx8 %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Filter out low quality reads using selected thresholds - these will change with experiment
utx8 <- subset(x = utx8, 
               subset= (nUMI >= 1000) & 
                 (nGene >= 500) & 
                 (log10GenesPerUMI > 0.83) & 
                 (mitoRatio < 0.20))

##################################   UTX13 #########
# Create metadata dataframe for wt8
metadata.utx13 <- utx13@meta.data

# Add cell IDs to metadata
metadata.utx13$cells <- rownames(metadata.utx13)

# Rename columns
metadata.utx13 <- metadata.utx13 %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata.utx13$sample <- NA
metadata.utx13$sample[which(str_detect(metadata.utx13$cells, "^utx596_"))] <- "utx596"
metadata.utx13$sample[which(str_detect(metadata.utx13$cells, "^utx600_"))] <- "utx600"
metadata.utx13$sample[which(str_detect(metadata.utx13$cells, "^utx13wk_"))] <- "utx13wk"

metadata.utx13$mice <- NA
metadata.utx13$mice <- "UTX"

metadata.utx13$week <- NA
metadata.utx13$week <- "Wk13_UTX"


# Add metadata back to Seurat object
utx13@meta.data <- metadata.utx13

metadata.utx13 %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata.utx13 %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata.utx13 %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata.utx13 %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata.utx13 %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata.utx13 %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata.utx13 %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Filter out low quality reads using selected thresholds - these will change with experiment
utx13 <- subset(x = utx13, 
                subset= (nUMI >= 1000) & 
                  (nGene >= 800) & 
                  (log10GenesPerUMI > 0.8) & 
                  (mitoRatio < 0.20))

# Gene-level filtering ####
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts.wt8 <- GetAssayData(object = wt8, slot = "counts")
counts.wt13 <- GetAssayData(object = wt13, slot = "counts")
counts.utx8 <- GetAssayData(object = utx8, slot = "counts")
counts.utx13 <- GetAssayData(object = utx13, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero.wt8 <- counts.wt8 > 0
nonzero.wt13 <- counts.wt13 > 0
nonzero.utx8 <- counts.utx8 > 0
nonzero.utx13 <- counts.utx13 > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes.wt8 <- Matrix::rowSums(nonzero.wt8) >= 10
keep_genes.wt13 <- Matrix::rowSums(nonzero.wt13) >= 10
keep_genes.utx8 <- Matrix::rowSums(nonzero.utx8) >= 10
keep_genes.utx13 <- Matrix::rowSums(nonzero.utx13) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts.wt8 <- counts.wt8[keep_genes.wt8, ]
filtered_counts.wt13 <- counts.wt13[keep_genes.wt13, ]
filtered_counts.utx8 <- counts.utx8[keep_genes.utx8, ]
filtered_counts.utx13 <- counts.utx13[keep_genes.utx13, ]

# Reassign to filtered Seurat object
wt8 <- CreateSeuratObject(filtered_counts.wt8, meta.data = wt8@meta.data)
wt13 <- CreateSeuratObject(filtered_counts.wt13, meta.data = wt13@meta.data)
utx8 <- CreateSeuratObject(filtered_counts.utx8, meta.data = utx8@meta.data)
utx13 <- CreateSeuratObject(filtered_counts.utx13, meta.data = utx13@meta.data)

# Cell cycle scoring ######### 
library(biomaRt)

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genes[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

m.s.genes <- convertHumanGeneList(cc.genes$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)
wt8 <- CellCycleScoring(wt8, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
wt13 <- CellCycleScoring(wt13, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
utx8 <- CellCycleScoring(utx8, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
utx13 <- CellCycleScoring(utx13, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

options(future.globals.maxSize = 20 * 1024^3)

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
DefaultAssay(wt8) <- "RNA"
DefaultAssay(wt13) <- "RNA"
DefaultAssay(utx8) <- "RNA"
DefaultAssay(utx13) <- "RNA"

wt8$CC.Difference <- wt8$S.Score - wt8$G2M.Score
wt13$CC.Difference <- wt13$S.Score - wt13$G2M.Score
utx8$CC.Difference <- utx8$S.Score - utx8$G2M.Score
utx13$CC.Difference <- utx13$S.Score - utx13$G2M.Score

wt8 <- NormalizeData(wt8)
wt8 <- FindVariableFeatures(wt8, selection.method = "vst", nfeatures = 2000)
wt13 <- NormalizeData(wt13)
wt13 <- FindVariableFeatures(wt13, selection.method = "vst", nfeatures = 2000)
utx8 <- NormalizeData(utx8)
utx8 <- FindVariableFeatures(utx8, selection.method = "vst", nfeatures = 2000)
utx13 <- NormalizeData(utx13)
utx13 <- FindVariableFeatures(utx13, selection.method = "vst", nfeatures = 2000)

split.list <- list(wt8,wt13,utx8,utx13)
ref.list <- split.list[c("WT","UTX")]


features <- SelectIntegrationFeatures(object.list = split.list)

split.list <- lapply(X = split.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

all.anchors <- FindIntegrationAnchors(object.list = split.list, reference = c(1,2), reduction = "rpca", 
                                      dims = 1:30)
all.integrated <- IntegrateData(anchorset = all.anchors, dims = 1:30)
DefaultAssay(all.integrated) <- "integrated"


all.integrated <- ScaleData(all.integrated, verbose = T, vars.to.regress=c("S.Score", "G2M.Score","percent.mt","nFeature_RNA","nCount_RNA"))
all.integrated <- RunPCA(all.integrated, npcs = 30, verbose = T)
all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:30)
all.integrated <- FindNeighbors(all.integrated, reduction = "pca", dims = 1:30)
all.integrated <- FindClusters(all.integrated, resolution = 0.35)

DimPlot(all.integrated, reduction = "umap", pt.size = 1,  label=T)
DimPlot(all.integrated, reduction = "umap", split.by = "mice", pt.size = 1,  label=T)
DimPlot(all.integrated, reduction = "umap", split.by = "week", pt.size = 1,  label=T)
DimPlot(all.integrated, reduction = "umap", group.by = "mice", pt.size = 1,  label=T)
Idents(all.integrated) <- "seurat_clusters"
Idents(all.integrated) <- factor(Idents(all.integrated), levels = c("Wk8_WT","Wk8_UTX","Wk13_WT","Wk13_UTX"))
head(all.integrated@meta.data)
DefaultAssay(all.integrated) <- "RNA"

#### Please modify based on your own resolution and cluster numbers
Idents(object = all.integrated) <- "integrated_snn_res.0.35"
current.cluster.ids <- c(0:16)
new.cluster.ids <- c("Cd4-1","Cd4-2",
                     "Cd8-1","Cd8-2",
                     "Cd4-3","B-1",
                     "Cd4-4","Cd4-5",
                     "DG","Cd4-6",
                     "B-2","Cd8-3",
                     "Cd8-4","Dying",
                     "Cd8-5","MC",
                     "B-3")

Idents(object = all.integrated) <- plyr::mapvalues(x = Idents(object = all.integrated), from = current.cluster.ids, to = new.cluster.ids)

## Store this new naming into a variable "$celltypes" in the metadata and define the ordering of the cell types
all.integrated@meta.data$celltypes <- Idents(object = all.integrated)
all.integrated@meta.data$celltypes <- factor(all.integrated@meta.data$celltypes,
                                             levels = c("Cd4-1","Cd4-2","Cd4-3","Cd4-4",
                                                        "Cd4-5","Cd4-6","Dying",
                                                        "Cd8-1","Cd8-2","Cd8-3","Cd8-4","Cd8-5",
                                                        "DG","B-1","B-2","B-3","MC"))
head(all.integrated@meta.data$celltypes)

Idents(all.integrated) <- "celltypes"
Idents(all.integrated) <- factor(Idents(all.integrated), levels = c("Cd4-1","Cd4-2","Cd4-3","Cd4-4",
                                                                    "Cd4-5","Cd4-6","Dying",
                                                                    "Cd8-1","Cd8-2","Cd8-3","Cd8-4","Cd8-5",
                                                                    "DG","B-1","B-2","B-3","MC"))


all.integrated <- readRDS("all.integrated.rds")

#markers
FeaturePlot(all.integrated, features = c("Cd3g","Ms4a1","Cd74"), label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)
FeaturePlot(all.integrated, features = c("H2-Oa","Btla", "Cd14"), order=T,label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)
FeaturePlot(all.integrated, features = c("Siglech","Trdc", "Ncr1"), order=T, label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)
FeaturePlot(all.integrated, features = c("Foxp3","Rorc", "Il2ra"), order=T, label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)
FeaturePlot(all.integrated, features = c("Cd8a", "Cd4","Ncr1"), order=T, label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)
FeaturePlot(all.integrated, features = c("Ncr1","Il17a"), order=T, label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)
FeaturePlot(all.integrated, features = c("Rora","Gata3"), order=T, label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)
FeaturePlot(all.integrated, features = c("Prdm1","Tnfrsf17"), order=T, label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)


FeaturePlot(all.integrated, features = c("Ccr9","Cxcr6","Itgb1"), order=T,label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)
FeaturePlot(all.integrated, features = c("Itga4","Itgb7"), order=T, label = FALSE, pt.size = 0.1, cols = c("grey", "red")) + DimPlot(all.integrated, reduction = "umap", label = TRUE, pt.size = 0.9)

#b: 5,10,16 - plasma =5?
#mac: 15 = myeloid
#no pdc
#8: gd t?
#4: foxp3
#cd4: 0,1,13,4,6,9,7 NOT 3
#cd8:2,3, 11,12,14
#nkt part of 8
#13: dying - cd4-7

Idents(all.integrated) <- "celltypes"
cd4.1 <- FindMarkers(all.integrated, ident.1 = "Cd4-1", logfc.threshold = 0.25, min.pct = 0)
cd4.2 <- FindMarkers(all.integrated, ident.1 = "Cd4-2", logfc.threshold = 0.25, min.pct = 0)
cd4.3 <- FindMarkers(all.integrated, ident.1 = "Cd4-3", logfc.threshold = 0.25, min.pct = 0)
cd4.4 <- FindMarkers(all.integrated, ident.1 = "Cd4-4", logfc.threshold = 0.25, min.pct = 0)
cd4.5 <- FindMarkers(all.integrated, ident.1 = "Cd4-5", logfc.threshold = 0.25, min.pct = 0)
cd4.6 <- FindMarkers(all.integrated, ident.1 = "Cd4-6", logfc.threshold = 0.25, min.pct = 0)
cd4.7 <- FindMarkers(all.integrated, ident.1 = "Cd4-7", logfc.threshold = 0.25, min.pct = 0)
cd8.1 <- FindMarkers(all.integrated, ident.1 = "Cd8-1", logfc.threshold = 0.25, min.pct = 0)
cd8.2 <- FindMarkers(all.integrated, ident.1 = "Cd8-2", logfc.threshold = 0.25, min.pct = 0)
cd8.3 <- FindMarkers(all.integrated, ident.1 = "Cd8-3", logfc.threshold = 0.25, min.pct = 0)
cd8.4 <- FindMarkers(all.integrated, ident.1 = "Cd8-4", logfc.threshold = 0.25, min.pct = 0)
cd8.5 <- FindMarkers(all.integrated, ident.1 = "Cd8-5", logfc.threshold = 0.25, min.pct = 0)
b.1 <- FindMarkers(all.integrated, ident.1 = "B-1", logfc.threshold = 0.25, min.pct = 0)
b.2 <- FindMarkers(all.integrated, ident.1 = "B-2", logfc.threshold = 0.25, min.pct = 0)
b.3 <- FindMarkers(all.integrated, ident.1 = "B-3", logfc.threshold = 0.25, min.pct = 0)
mc <- FindMarkers(all.integrated, ident.1 = "MC", logfc.threshold = 0.25, min.pct = 0)
dg <- FindMarkers(all.integrated, ident.1 = "DG", logfc.threshold = 0.25, min.pct = 0)

write.csv(all.1, "DEG_ydens_1.csv")
write.csv(all.2, "DEG_ydens_2.csv")
write.csv(all.3, "DEG_ydens_3.csv")
write.csv(all.4, "DEG_ydens_4.csv")
write.csv(all.5, "DEG_ydens_5.csv")
write.csv(all.6, "DEG_ydens_6.csv")
write.csv(all.7, "DEG_ydens_7.csv")
write.csv(all.8, "DEG_ydens_8.csv")
write.csv(all.9, "DEG_ydens_9.csv")
write.csv(all.10, "DEG_ydens_10.csv")
write.csv(all.11, "DEG_ydens_11.csv")
write.csv(all.12, "DEG_ydens_12.csv")
write.csv(all.13, "DEG_ydens_13.csv")
write.csv(all.14, "DEG_ydens_14.csv")
write.csv(all.15, "DEG_ydens_15.csv")
write.csv(all.16, "DEG_ydens_16.csv")

library(msigdbr)
library(fgsea)

#msigdb
msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

ranks <- cd8.5$avg_logFC
names(ranks) <- rownames(cd8.5)

# allPathways is predefined list of signatures
fgseaRes <- fgsea(pathways = fgsea_sets, 
                  # scoreType = "pos",
                  stats = ranks,
                  minSize=10,
                  maxSize=500,
                  eps=0
)

head(fgseaRes[order(padj), ])

ggplot(fgseaRes %>% filter(padj < 0.05) %>% head(n=10), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways in Cd8-5") + 
  theme_minimal()

heatmap

Idents(all.integrated) <- "celltypes"
Idents(all.integrated) <- "mice"
Idents(all.integrated) <- "week"

markers.tnf=c("Tnf","Tnfrsf1a","Tnfrsf1b")
VlnPlot(all.integrated, markers.tnf, pt.size = 0)

markers.ifng=c("Ifng","Ifngr1","Ifngr2")
VlnPlot(all.integrated, markers.ifng, pt.size = 0)

markers.ifna=c("Ifna","Ifnar1","Ifnar2")
VlnPlot(all.integrated, markers.ifna, pt.size = 0)

markers.chemo1=c("Ccr9","Cxcr6","Itgb1")
VlnPlot(all.integrated, markers.chemo1pt.size = 0)

markers.ccr9=c("Ccr9")
VlnPlot(all.integrated, markers.ccr9, split.by = "week", pt.size = 0)
markers.cxcr6=c("Cxcr6")
VlnPlot(all.integrated, markers.cxcr6, split.by = "week", pt.size = 0)
markers.itgb1=c("Itgb1")
VlnPlot(all.integrated, markers.itgb1, split.by = "week", pt.size = 0)
markers.itgb1=c("Cx3cr1")
VlnPlot(all.integrated, markers.itgb1, split.by = "week", pt.size = 0)
markers.itgb1=c("Ccr4")
VlnPlot(all.integrated, markers.itgb1, split.by = "week", pt.size = 0)


markers.chemo2=c("Itga4","Itgb7")
VlnPlot(all.integrated, markers.chemo2, pt.size = 0)

head(all.merged@meta.data)
saveRDS(all.merged, "merged_lisa.rds")


# How does cluster membership vary by genotype?
Idents(all.integrated) <- "celltypes"
clustercount <- table(Idents(all.integrated), all.integrated$week)
write.csv(clustercount, "lisa_weeks-count-by-celltypes.csv")

all.prop.m <- read.csv("lisanumbers_m.csv")
all.prop.m$celltypes <- factor(all.prop.m$celltypes,levels = c("Cd4-1", "Cd4-2", "Cd4-3", "Cd4-4","Cd4-5","Cd4-6",
                                                         "Cd8-1","Cd8-2","Cd8-3","Cd8-4","Cd8-5",
                                                         "DG",
                                                         "B-1","B-2","B-3",
                                                         "MC"))
all.prop.m$Mice <- factor(all.prop.m$Mice,levels = c("WT Week 8","UTX Week 8","WT Week 13","UTX Week 13"))

head(all.prop.m)
all.prop_m <- ggplot(all.prop.m, aes(x = Mice, y = prop, fill = celltypes)) +
  theme_bw(base_size = 10) +
  geom_col(position = "stack", width = 0.7) +
  xlab("Mice") +
  ylab("Proportion of immune cells") +
  theme(legend.text = element_text(size = 9), 
        legend.title = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(angle = 90, size = 11), 
        axis.title = element_text(face = "bold", size = 11))

all.c_m = factor(c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))
all.c.m <- colorRampPalette(all.c_m)

all.prop_m + labs(stack = "Mice") + scale_fill_manual(values = all.c.m(18))


all.prop.w <- read.csv("lisanumbers_w.csv")
all.prop.w$Weeks <- factor(all.prop.w$Weeks,levels = c("Week 8","Week 13"))
all.prop.w$celltypes <- factor(all.prop.w$celltypes,levels = c("Cd4-1", "Cd4-2", "Cd4-3", "Cd4-4","Cd4-5","Cd4-6",
                                                         "Cd8-1","Cd8-2","Cd8-3","Cd8-4","Cd8-5",
                                                         "DG",
                                                         "B-1","B-2","B-3",
                                                         "MC"))

head(all.prop.w)

library(wesanderson)
all.prop_color <- brewer.pal(n =12, name = "Paired")

b <- ggplot(all.prop.m, aes(x = celltypes, y = prop)) +
  geom_bar(stat = "identity", position = position_dodge(), aes(fill=Mice)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) + labs(title="Proportion of immune cells", x ="Clusters", y = "Proportion")


all.prop_w <- ggplot(all.prop.w, aes(x = Weeks, y = prop, fill = celltypes)) +
  theme_bw(base_size = 10) +
  geom_col(position = "stack", width = 0.7) +
  xlab("Mice") +
  ylab("Proportion of immune cells") +
  theme(legend.text = element_text(size = 9), 
        legend.title = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(angle = 90, size = 11), 
        axis.title = element_text(face = "bold", size = 11))


all.prop_w + labs(stack = "Weeks") + scale_fill_manual(values = all.c.m(18))

library(ggpubr)

# Bar chart separated by cell type
ggplot(all.prop.w, aes(x = celltypes, y = prop, fill = Weeks))+
geom_bar(stat = "identity", position = position_dodge())+
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
scale_fill_manual(values = all.c.m(18))+ labs(title="Proportion of immune cells", x ="Clusters", y = "Proportion")+
stat_compare_means(aes(group = Weeks), label = "p.signif", label.y = c(0.3,0.23,0.1,0.06,0.06,0.04,0.15,0.11,0.03,0.03,0.03,0.03,0.1,0.04,0.02,0.02))

library(rstatix)
# Transform `dose` into factor variable
df <- all.prop.m
df$Mice <- as.factor(df$Mice)


stat.test <- df %>%
  group_by(celltypes) %>%
  t_test(prop ~ Mice)
stat.test
stat.test1 <- stat.test %>% add_xy_position(x = "celltypes")

# Bar chart separated by cell type
ggplot(all.prop.m, aes(x = celltypes, y = prop))+
  geom_bar(stat = "identity", position = position_dodge(), aes(fill = Mice))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) + 
  labs(title="Proportion of immune cells", x ="Clusters", y = "Proportion") +
  stat_pvalue_manual(
  stat.test1, label = "p.adj.signif", tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.4)))
