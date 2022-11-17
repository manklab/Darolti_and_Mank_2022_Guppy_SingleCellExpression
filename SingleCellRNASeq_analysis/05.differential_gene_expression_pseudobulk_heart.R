
heart_umap_rename_clusters_dbrmv_PB <- heart_umap_rename_clusters_dbrmv
heart_umap_rename_clusters_dbrmv_PB$CellType <- "All"

counts_heart_PB <- heart_umap_rename_clusters_dbrmv_PB@assays$RNA@counts
metadata_heart_PB <- heart_umap_rename_clusters_dbrmv_PB@meta.data
metadata_heart_PB$CellType <- "All"
sce_heart_PB <- SingleCellExperiment(assays=list(counts=counts_heart_PB), colData=metadata_heart_PB)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_heart_PB)$sample_id <- as.factor(colData(sce_heart_PB)$sample)
colData(sce_heart_PB)$sex_id <- as.factor(colData(sce_heart_PB)$sex)
colData(sce_heart_PB)$cluster_id <- as.factor(colData(sce_heart_PB)$CellType)
kids_heart_PB <- purrr::set_names(levels(sce_heart_PB$cluster_id))
nk_heart_PB <- length(kids_heart_PB)
sids_heart_PB <- purrr::set_names(levels(sce_heart_PB$sample_id))
ns_heart_PB <- length(sids_heart_PB)

# Turn named vector into a numeric vector
n_cells_heart_PB <- as.numeric(table(sce_heart_PB$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_heart_PB <- match(sids_heart_PB, sce_heart_PB$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_heart_PB <- data.frame(colData(sce_heart_PB)[m_heart_PB, ], n_cells_heart_PB, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_heart_PB <- sce_heart_PB[rowSums(counts(sce_heart_PB)>1)>=10,]

# Aggregate counts per sample_id and cluster_id
groups_heart_PB <- colData(sce_heart_PB)[, c("cluster_id", "sample_id")]
pb_heart_PB <- aggregate.Matrix(t(counts(sce_heart_PB)), groupings=groups_heart_PB, fun="sum")

# Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_heart_PB <- sapply(stringr::str_split(rownames(pb_heart_PB), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_heart_PB <- split.data.frame(pb_heart_PB, factor(splitf_heart_PB)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# Get sample names
get_sample_ids_heart_PB <- function(x){
	pb_heart_PB[[x]] %>% colnames()
}
de_samples_heart_PB <- map(1:length(kids_heart_PB), get_sample_ids_heart_PB) %>% unlist()

# Get cluster IDs (in this case there is just one as we are merging counts across all cells)
samples_list_heart_PB <- map(1:length(kids_heart_PB), get_sample_ids_heart_PB)
get_cluster_ids_heart_PB <- function(x){
	rep(names(pb_heart_PB)[x], each=length(samples_list_heart_PB[[x]]))
}
de_cluster_ids_heart_PB <- map(1:length(kids_heart_PB), get_cluster_ids_heart_PB) %>% unlist()

# Create dataframe
gg_df_heart_PB <- data.frame(cluster_id=de_cluster_ids_heart_PB, sample_id=de_samples_heart_PB)
gg_df_heart_PB <- left_join(gg_df_heart_PB, ei_heart_PB[, c("sample_id", "sex_id", "n_cells_heart_PB")])
metadata_heart_PB <- gg_df_heart_PB %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_heart_PB)

metadata_heart_PB$cluster_id <- as.factor(metadata_heart_PB$cluster_id)

clusters_heart_PB <- levels(metadata_heart_PB$cluster_id)

cluster_metadata_heart_PB <- metadata_heart_PB[which(metadata_heart_PB$cluster_id == clusters_heart_PB[1]), ]
rownames(cluster_metadata_heart_PB) <- cluster_metadata_heart_PB$sample_id
counts_heart_PB <- pb_heart_PB[[clusters_heart_PB[1]]]
cluster_counts_heart_PB <- data.frame(counts_heart_PB[, which(colnames(counts_heart_PB) %in% rownames(cluster_metadata_heart_PB))])
all(rownames(cluster_metadata_heart_PB) == colnames(cluster_counts_heart_PB))    
		[1] TRUE

# Create DESeq2 object
dds_heart_PB <- DESeqDataSetFromMatrix(cluster_counts_heart_PB, 
				colData = cluster_metadata_heart_PB, 
                design = ~ sex_id)
raw_counts_PB <- counts(dds_heart_PB)
write.csv(raw_counts_PB, paste0("~/bulk/dds_raw_counts_PB_heart.csv"), quote = FALSE, row.names = T)
cpm_PB <- calculateCPM(dds_heart_PB)
write.csv(cpm_PB, paste0("~/bulk/dds_cpm_PB_heart.csv"), quote = FALSE, row.names = T)

# Transform counts for data visualization
rld_heart_PB <- rlog(dds_heart_PB, blind=TRUE)

# Run DESeq2 differential expression analysis
dds_deseq_heart_PB <- DESeq(dds_heart_PB)

# Output results of Wald test for contrast A vs B
contrast_heart_PB <- c("sex_id", levels(cluster_metadata_heart_PB$sex_id)[2], levels(cluster_metadata_heart_PB$sex_id)[1])
        [1] "sex_id" "Male"   "Female"
# resultsNames(dds)
res_heart_PB <- results(dds_deseq_heart_PB, 
                contrast = contrast_heart_PB,
                alpha = 0.05)
                      
# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1
        
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_heart_PB <- res_heart_PB %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
        
# Subset the significant results
sig_res_heart_PB <- dplyr::filter(res_tbl_heart_PB, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_heart_PB, paste0("~/bulk/deseq_sig_genes_heart.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_heart <- dplyr::filter(sig_res_heart_PB, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_heart <- dplyr::filter(sig_res_heart_PB, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_heart, paste0("~/bulk/deseq_sig_genes_fcpass_malebiased_heart.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_heart, paste0("~/bulk/deseq_sig_genes_fcpass_femalebiased_heart.csv"), quote = FALSE, row.names = FALSE)

