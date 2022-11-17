
skin_umap_rename_clusters_dbrmv_PB <- skin_umap_rename_clusters_dbrmv
skin_umap_rename_clusters_dbrmv_PB$CellType <- "All"

counts_skin_PB <- skin_umap_rename_clusters_dbrmv_PB@assays$RNA@counts
metadata_skin_PB <- skin_umap_rename_clusters_dbrmv_PB@meta.data
metadata_skin_PB$CellType <- "All"
sce_skin_PB <- SingleCellExperiment(assays=list(counts=counts_skin_PB), colData=metadata_skin_PB)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_skin_PB)$sample_id <- as.factor(colData(sce_skin_PB)$sample)
colData(sce_skin_PB)$sex_id <- as.factor(colData(sce_skin_PB)$sex)
colData(sce_skin_PB)$cluster_id <- as.factor(colData(sce_skin_PB)$CellType)
kids_skin_PB <- purrr::set_names(levels(sce_skin_PB$cluster_id))
nk_skin_PB <- length(kids_skin_PB)
sids_skin_PB <- purrr::set_names(levels(sce_skin_PB$sample_id))
ns_skin_PB <- length(sids_skin_PB)

# Turn named vector into a numeric vector
n_cells_skin_PB <- as.numeric(table(sce_skin_PB$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_skin_PB <- match(sids_skin_PB, sce_skin_PB$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_skin_PB <- data.frame(colData(sce_skin_PB)[m_skin_PB, ], n_cells_skin_PB, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_skin_PB <- sce_skin_PB[rowSums(counts(sce_skin_PB)>1)>=10,]

# Aggregate counts per sample_id and cluster_id
groups_skin_PB <- colData(sce_skin_PB)[, c("cluster_id", "sample_id")]
pb_skin_PB <- aggregate.Matrix(t(counts(sce_skin_PB)), groupings=groups_skin_PB, fun="sum")

# Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_skin_PB <- sapply(stringr::str_split(rownames(pb_skin_PB), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_skin_PB <- split.data.frame(pb_skin_PB, factor(splitf_skin_PB)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# Get sample names
get_sample_ids_skin_PB <- function(x){
	pb_skin_PB[[x]] %>% colnames()
}
de_samples_skin_PB <- map(1:length(kids_skin_PB), get_sample_ids_skin_PB) %>% unlist()

# Get cluster IDs (in this case there is just one as we are merging counts across all cells)
samples_list_skin_PB <- map(1:length(kids_skin_PB), get_sample_ids_skin_PB)
get_cluster_ids_skin_PB <- function(x){
	rep(names(pb_skin_PB)[x], each=length(samples_list_skin_PB[[x]]))
}
de_cluster_ids_skin_PB <- map(1:length(kids_skin_PB), get_cluster_ids_skin_PB) %>% unlist()

# Create dataframe
gg_df_skin_PB <- data.frame(cluster_id=de_cluster_ids_skin_PB, sample_id=de_samples_skin_PB)
gg_df_skin_PB <- left_join(gg_df_skin_PB, ei_skin_PB[, c("sample_id", "sex_id", "n_cells_skin_PB")])
metadata_skin_PB <- gg_df_skin_PB %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_skin_PB)

metadata_skin_PB$cluster_id <- as.factor(metadata_skin_PB$cluster_id)

clusters_skin_PB <- levels(metadata_skin_PB$cluster_id)

cluster_metadata_skin_PB <- metadata_skin_PB[which(metadata_skin_PB$cluster_id == clusters_skin_PB[1]), ]
rownames(cluster_metadata_skin_PB) <- cluster_metadata_skin_PB$sample_id
counts_skin_PB <- pb_skin_PB[[clusters_skin_PB[1]]]
cluster_counts_skin_PB <- data.frame(counts_skin_PB[, which(colnames(counts_skin_PB) %in% rownames(cluster_metadata_skin_PB))])
all(rownames(cluster_metadata_skin_PB) == colnames(cluster_counts_skin_PB))    
		[1] TRUE

# Create DESeq2 object
dds_skin_PB <- DESeqDataSetFromMatrix(cluster_counts_skin_PB, 
				colData = cluster_metadata_skin_PB, 
                design = ~ sex_id)
raw_counts_PB <- counts(dds_skin_PB)
write.csv(raw_counts_PB, paste0("~/bulk/dds_raw_counts_PB_skin.csv"), quote = FALSE, row.names = T)
cpm_PB <- calculateCPM(dds_skin_PB)
write.csv(cpm_PB, paste0("~/bulk/dds_cpm_PB_skin.csv"), quote = FALSE, row.names = T)

# Transform counts for data visualization
rld_skin_PB <- rlog(dds_skin_PB, blind=TRUE)

# Run DESeq2 differential expression analysis
dds_deseq_skin_PB <- DESeq(dds_skin_PB)

# Output results of Wald test for contrast A vs B
contrast_skin_PB <- c("sex_id", levels(cluster_metadata_skin_PB$sex_id)[2], levels(cluster_metadata_skin_PB$sex_id)[1])
        [1] "sex_id" "Male"   "Female"
# resultsNames(dds)
res_skin_PB <- results(dds_deseq_skin_PB, 
                contrast = contrast_skin_PB,
                alpha = 0.05)
                      
# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1
        
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_skin_PB <- res_skin_PB %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
        
# Subset the significant results
sig_res_skin_PB <- dplyr::filter(res_tbl_skin_PB, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_skin_PB, paste0("~/bulk/deseq_sig_genes_skin.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_skin <- dplyr::filter(sig_res_skin_PB, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_skin <- dplyr::filter(sig_res_skin_PB, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_skin, paste0("~/bulk/deseq_sig_genes_fcpass_malebiased_skin.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_skin, paste0("~/bulk/deseq_sig_genes_fcpass_femalebiased_skin.csv"), quote = FALSE, row.names = FALSE)

