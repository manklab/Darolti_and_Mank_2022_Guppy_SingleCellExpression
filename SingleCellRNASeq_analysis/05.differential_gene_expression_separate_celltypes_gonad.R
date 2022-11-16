library(SingleCellExperiment)
library(Seurat)
library(Matrix.utils)
library(Matrix)
library(dplyr)
library(magrittr)
library(purrr)
library(stringr)
library(DESeq2)
library(pheatmap)
library(apeglm)
library(tibble)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scran)
library(scuttle)

load("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/03.marker_identification/gonad_filtered_dbrmv_reevaluated/gonad_umap_3_rename_clusters.RData")

head(gonad_umap_3_rename_clusters)

counts_gonad <- gonad_umap_3_rename_clusters@assays$RNA@counts
metadata_gonad <- gonad_umap_3_rename_clusters@meta.data
metadata_gonad$CellType <- factor(gonad_umap_3_rename_clusters@active.ident)

gonad_subset <- subset(gonad_umap_3_rename_clusters, idents="Germ", invert=T)
counts_gonad <- gonad_subset@assays$RNA@counts
metadata_gonad <- gonad_subset@meta.data
metadata_gonad$CellType <- factor(gonad_subset@active.ident)

sce_gonad <- SingleCellExperiment(assays=list(counts=counts_gonad), colData=metadata_gonad)
sce_gonad 
class: SingleCellExperiment 
dim: 21793 26695 
metadata(0):
assays(1): counts
rownames(21793): ENSPREG00000006268 ENSPREG00000006273 ... ENSPREG00000011737 BDKRB2.2
rowData names(0):
colnames(26695): F1_Lib11_AAACCCAAGGTACAGC-1 F1_Lib11_AAACCCAAGTCTGGAG-1 ...
  M6_Lib37_TTTGTTGTCGCACGAC-1 M6_Lib37_TTTGTTGTCTATACTC-1
colData names(17): seq_folder nUMI ... seurat_clusters CellType
reducedDimNames(0):
altExpNames(0):

assays(sce_gonad)
		List of length 1
		names(1): counts
dim(counts(sce_gonad))
[1] 21793 26695
counts(sce_gonad)[1:2, 1:2]
		2 x 2 sparse Matrix of class "dgCMatrix"
		                   F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1
		ENSPREG00000006268                          .                          .
		arl13b                                      .                          .
dim(colData(sce_gonad))
[1] 26695    17

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_gonad)$sample_id <- as.factor(colData(sce_gonad)$sample)
colData(sce_gonad)$sex_id <- as.factor(colData(sce_gonad)$sex)
colData(sce_gonad)$cluster_id <- as.factor(colData(sce_gonad)$CellType)

kids_gonad <- purrr::set_names(levels(sce_gonad$cluster_id))
kids_gonad
        Spermatocyte          Spermatozoa            Granulosa           Macrophage 
      "Spermatocyte"        "Spermatozoa"          "Granulosa"         "Macrophage" 
           Endocrine   Mitochondrial-rich   Embryonic Yolk Sac 
         "Endocrine" "Mitochondrial-rich" "Embryonic Yolk Sac"

nk_gonad <- length(kids_gonad)
nk_gonad
		[1] 7

sids_gonad <- purrr::set_names(levels(sce_gonad$sample_id))
sids_gonad
		   GonadF1   GonadF2   GonadF3   GonadM4   GonadM5   GonadM6 
		 "GonadF1" "GonadF2" "GonadF3" "GonadM4" "GonadM5" "GonadM6" 
ns_gonad <- length(sids_gonad)
ns_gonad
		[1] 6

# Number of cells per sample
table(sce_gonad$sample_id)
GonadF1 GonadF2 GonadF3 GonadM1 GonadM2 GonadM3 
   2572    1635    1191    7973    6541    6783 
  
# Turn named vector into a numeric vector
n_cells_gonad <- as.numeric(table(sce_gonad$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_gonad <- match(sids_gonad, sce_gonad$sample_id)
m_gonad
[1]     1  2573  4208  5399 13372 19913

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_gonad <- data.frame(colData(sce_gonad)[m_gonad, ], n_cells_gonad, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
dim(sce_gonad)
[1] 21793 26695
sce_gonad <- sce_gonad[rowSums(counts(sce_gonad)>0)>0,]
dim(sce_gonad)	
[1] 20912 26695
sce_gonad <- sce_gonad[rowSums(counts(sce_gonad)>1)>=10,]
dim(sce_gonad)
[1] 12994 26695

# Aggregate counts per sample_id and cluster_id
groups_gonad <- colData(sce_gonad)[, c("cluster_id", "sample_id")]
pb_gonad <- aggregate.Matrix(t(counts(sce_gonad)), groupings=groups_gonad, fun="sum")
class(pb_gonad)
		[1] "dgCMatrix"
		attr(,"package")
		[1] "Matrix"
dim(pb_gonad)
[1]    41 12994

pb_gonad[1:6, 1:6]
		6 x 6 sparse Matrix of class "dgCMatrix"
                           arl13b zgc:152904 gart n6amt1 timmdc1 VANGL1
Embryonic Yolk Sac_GonadF1      .          4    .      2       .      1
Embryonic Yolk Sac_GonadF2      1          1    2      4       .      .
Embryonic Yolk Sac_GonadF3      1          .    .      .       .      .
Embryonic Yolk Sac_GonadM4      .          1    .      1       2      .
Embryonic Yolk Sac_GonadM5      1          .    .      .       .      .
Embryonic Yolk Sac_GonadM6      .          .    .      .       .      .


## Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_gonad <- sapply(stringr::str_split(rownames(pb_gonad), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_gonad <- split.data.frame(pb_gonad, factor(splitf_gonad)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
class(pb_gonad)
		[1] "list"
str(pb_gonad)

# Number of cells for each cell type and sample
options(width=100)
table(sce_gonad$cluster_id, sce_gonad$sample_id)
                     GonadF1 GonadF2 GonadF3 GonadM4 GonadM5 GonadM6
  Spermatocyte            27      55      29    7324    4838    5633
  Spermatozoa              2       2       0     324    1402    1024
  Granulosa              648     455     366      89      11      14
  Macrophage             881     268     260      18      38      13
  Endocrine              566     518     305     147     235      81
  Mitochondrial-rich     357     234     154      11       2       4
  Embryonic Yolk Sac      91     103      77      60      15      14


get_sample_ids_gonad <- function(x){
	pb_gonad[[x]] %>% colnames()
}
de_samples_gonad <- map(1:length(kids_gonad), get_sample_ids_gonad) %>% unlist()

samples_list_gonad <- map(1:length(kids_gonad), get_sample_ids_gonad)

get_cluster_ids_gonad <- function(x){
	rep(names(pb_gonad)[x], each=length(samples_list_gonad[[x]]))
}

de_cluster_ids_gonad <- map(1:length(kids_gonad), get_cluster_ids_gonad) %>% unlist()

gg_df_gonad <- data.frame(cluster_id=de_cluster_ids_gonad, sample_id=de_samples_gonad)
gg_df_gonad <- left_join(gg_df_gonad, ei_gonad[, c("sample_id", "sex_id", "n_cells_gonad")])
metadata_gonad <- gg_df_gonad %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_gonad)

metadata_gonad$cluster_id <- as.factor(metadata_gonad$cluster_id)

clusters_gonad <- levels(metadata_gonad$cluster_id)

get_dds_resultsAvsB <- function(x, A, B){
        cluster_metadata_gonad <- metadata_gonad[which(metadata_gonad$cluster_id == clusters_gonad[x]), ]
        rownames(cluster_metadata_gonad) <- cluster_metadata_gonad$sample_id
        counts_gonad <- pb_gonad[[clusters_gonad[x]]]
        cluster_counts_gonad <- data.frame(counts_gonad[, which(colnames(counts_gonad) %in% rownames(cluster_metadata_gonad))])
        
        all(rownames(cluster_metadata_gonad) == colnames(cluster_counts_gonad))        
        
        cpm_gonad <- calculateCPM(cluster_counts_gonad)
        write.csv(cpm_gonad, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_gonad[x], "_dds_cpm_test.csv"), quote = FALSE, row.names = T)
        

        
       		
}      

map(1:length(clusters_gonad), get_dds_resultsAvsB, A = 2, B = 1)


str("Create DESeq2 object")
        dds_gonad <- DESeqDataSetFromMatrix(cluster_counts_gonad, 
                                      colData = cluster_metadata_gonad, 
                                      design = ~ sex_id)
        raw_counts <- counts(dds_gonad)
        write.csv(raw_counts, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/separate_clusters/",clusters_gonad[x],"_dds_raw_counts.csv"), quote = FALSE, row.names = T)
        
        fpm_gonad <- fpm(dds_gonad)
        write.csv(fpm_gonad, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_gonad[x], "_dds_fpm.csv"), quote = FALSE, row.names = T)
        
        cpm_gonad <- calculateCPM(dds_gonad)
        write.csv(cpm_gonad, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_gonad[x], "_dds_cpm.csv"), quote = FALSE, row.names = T)
        




