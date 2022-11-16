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
library(scuttle)
library(ashr)

load("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/04.remove_doublets/heart_filtered/heart_umap_3_rename_clusters_dbrmv.RData")

heart_umap_3_rename_clusters_dbrmv_PB <- heart_umap_3_rename_clusters_dbrmv
heart_umap_3_rename_clusters_dbrmv_PB$CellType <- "All"

head(heart_umap_3_rename_clusters_dbrmv_PB)

counts_heart_PB <- heart_umap_3_rename_clusters_dbrmv_PB@assays$RNA@counts
metadata_heart_PB <- heart_umap_3_rename_clusters_dbrmv_PB@meta.data
metadata_heart_PB$CellType <- "All"

sce_heart_PB <- SingleCellExperiment(assays=list(counts=counts_heart_PB), colData=metadata_heart_PB)
sce_heart_PB
class: SingleCellExperiment  
dim: 21975 28923 
metadata(0):
assays(1): counts
rownames(21975): ENSPREG00000006268 arl13b ... guca1d ENSPREG00000006408
rowData names(0):
colnames(28923): F1_Lib6_AAACCCACAACCTATG-1 F1_Lib6_AAACCCACAACGGCTC-1 ...
  M3_Lib27_TTTGGTTGTATATGGA-1 M3_Lib27_TTTGGTTGTGTTACTG-1
colData names(27): seq_folder nUMI ... seurat_clusters CellType
reducedDimNames(0):
altExpNames(0):

assays(sce_heart_PB)
		List of length 1
		names(1): counts
dim(counts(sce_heart_PB))
[1] 21975 28923
counts(sce_heart_PB)[1:2, 1:2]
		2 x 2 sparse Matrix of class "dgCMatrix"
		                   F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1
		ENSPREG00000006268                          .                          .
		arl13b                                      .                          .
dim(colData(sce_heart_PB))
[1] 28923    27

# Determine the number of clusters and the cluster names
colData(sce_heart_PB)$sample_id <- as.factor(colData(sce_heart_PB)$sample)
colData(sce_heart_PB)$sex_id <- as.factor(colData(sce_heart_PB)$sex)
colData(sce_heart_PB)$cluster_id <- as.factor(colData(sce_heart_PB)$CellType)
kids_heart_PB <- purrr::set_names(levels(sce_heart_PB$cluster_id))
kids_heart_PB
		  All 
		"All" 
nk_heart_PB <- length(kids_heart_PB)
nk_heart_PB
		[1] 1

sids_heart_PB <- purrr::set_names(levels(sce_heart_PB$sample_id))
sids_heart_PB
		  HeartF1   HeartF2   HeartF3   HeartM1   HeartM2   HeartM3 
		"HeartF1" "HeartF2" "HeartF3" "HeartM1" "HeartM2" "HeartM3" 
ns_heart_PB <- length(sids_heart_PB)
ns_heart_PB
		[1] 6

# Determine no. cells per sample
table(sce_heart_PB$sample_id)
		HeartF1 HeartF2 HeartF3 HeartM1 HeartM2 HeartM3 
   6735    8966    6150    4615    1163    1294

# Turn named vector into a numeric vector
n_cells_heart_PB <- as.numeric(table(sce_heart_PB$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_heart_PB <- match(sids_heart_PB, sce_heart_PB$sample_id)
m_heart_PB
[1]     1  6736 15702 21852 26467 27630

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_heart_PB <- data.frame(colData(sce_heart_PB)[m_heart_PB, ], n_cells_heart_PB, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
dim(sce_heart_PB)
[1] 21975 28923
sce_heart_PB <- sce_heart_PB[rowSums(counts(sce_heart_PB)>1)>=10,]
dim(sce_heart_PB)
[1]  9860 28923

# Aggregate counts per sample_id and cluster_id
groups_heart_PB <- colData(sce_heart_PB)[, c("cluster_id", "sample_id")]
pb_heart_PB <- aggregate.Matrix(t(counts(sce_heart_PB)), groupings=groups_heart_PB, fun="sum")
class(pb_heart_PB)
		[1] "dgCMatrix"
		attr(,"package")
		[1] "Matrix"
dim(pb_heart_PB)
[1]     6 11478

pb_heart_PB[1:6, 1:6]
6 x 6 sparse Matrix of class "dgCMatrix"
            zgc:152904 gart VANGL1 sap18 rtraf gpr137c
All_HeartF1        217   43     80   887   132     164
All_HeartF2        654   40    105   911   155     192
All_HeartF3       1126   44    158   442    88     278
All_HeartM1        208   16     48   276    63      91
All_HeartM2        253    8     53   116    24      35
All_HeartM3        387    8     58   154    21     106


## Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_heart_PB <- sapply(stringr::str_split(rownames(pb_heart_PB), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_heart_PB <- split.data.frame(pb_heart_PB, factor(splitf_heart_PB)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb_heart_PB)
		[1] "list"
str(pb_heart_PB)

options(width=100)
table(sce_heart_PB$cluster_id, sce_heart_PB$sample_id)
		      HeartF1 HeartF2 HeartF3 HeartM1 HeartM2 HeartM3
  All    6735    8966    6150    4615    1163    1294
  
get_sample_ids_heart_PB <- function(x){
	pb_heart_PB[[x]] %>% colnames()
}
de_samples_heart_PB <- map(1:length(kids_heart_PB), get_sample_ids_heart_PB) %>% unlist()

samples_list_heart_PB <- map(1:length(kids_heart_PB), get_sample_ids_heart_PB)

get_cluster_ids_heart_PB <- function(x){
	rep(names(pb_heart_PB)[x], each=length(samples_list_heart_PB[[x]]))
}

de_cluster_ids_heart_PB <- map(1:length(kids_heart_PB), get_cluster_ids_heart_PB) %>% unlist()

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
write.csv(raw_counts_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/dds_raw_counts_PB.csv"), quote = FALSE, row.names = T)

#-----------
fpm_PB <- fpm(dds_heart_PB)
head(fpm_PB)
             HeartF1   HeartF2   HeartF3   HeartM1   HeartM2    HeartM3
zgc:152904 12.788239 37.472967 98.336041 33.946243 77.182608 103.985758
gart        2.534075  2.291925  3.842616  2.611249  2.440557   2.149576
VANGL1      4.714558  6.016302 13.798485  7.833748 16.168689  15.584429
sap18      52.272664 52.198583 38.600826 45.044053 35.388073  41.379346
rtraf       7.779021  8.881208  7.685232 10.281795  7.321670   5.642638
gpr137c     9.664844 11.001238 24.278348 14.851481 10.677436  28.481887

write.csv(fpm_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/dds_fpm_PB.csv"), quote = FALSE, row.names = T)

cpm_PB <- calculateCPM(dds_heart_PB)
head(cpm_PB)
             HeartF1   HeartF2   HeartF3   HeartM1   HeartM2   HeartM3
zgc:152904 15.151939 35.630845 98.938772 31.126732 93.091609 89.997342
gart        3.002458  2.179257  3.866169  2.394364  2.943608  1.860410
VANGL1      5.585968  5.720549 13.883060  7.183092 19.501404 13.487974
sap18      61.934425 49.632569 38.837422 41.302779 42.682319 35.812896
rtraf       9.216848  8.444619  7.732337  9.427808  8.830825  4.883577
gpr137c    11.451235 10.460432 24.427157 13.617945 12.878286 24.650435

write.csv(cpm_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/dds_cpm_PB.csv"), quote = FALSE, row.names = T)
#--------------


# Transform counts for data visualization
rld_heart_PB <- rlog(dds_heart_PB, blind=TRUE)
colData(rld_heart_PB)
DataFrame with 6 rows and 5 columns
        cluster_id   sample_id   sex_id n_cells_heart_PB sizeFactor
          <factor> <character> <factor>        <numeric>  <numeric>
HeartF1        All     HeartF1   Female             6735   2.162186
HeartF2        All     HeartF2   Female             8966   2.223841
HeartF3        All     HeartF3   Female             6150   1.459048
HeartM1        All     HeartM1   Male               4615   0.780756
HeartM2        All     HeartM2   Male               1163   0.417681
HeartM3        All     HeartM3   Male               1294   0.474221
        
# Plot PCA
DESeq2::plotPCA(rld_heart_PB, ntop=100,intgroup = "sex_id")
ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_PCAplot.png"))

        
        DESeq2::plotPCA(rld_heart_PB, ntop=100, intgroup = "n_cells_heart_PB")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_PCAplot_cellcount.png"))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat_heart_PB <- assay(rld_heart_PB)
head(rld_mat_heart_PB)
            HeartF1  HeartF2  HeartF3  HeartM1  HeartM2  HeartM3
zgc:152904 7.547091 8.369423 9.253767 8.288064 9.018492 9.305320
gart       4.346327 4.294971 4.580667 4.362242 4.331355 4.274497
VANGL1     5.744411 5.905791 6.540243 6.098664 6.662446 6.633398
sap18      8.582973 8.581708 8.314742 8.449653 8.243240 8.375517
rtraf      5.939212 6.038410 5.930405 6.146536 5.898020 5.728457
gpr137c    6.563636 6.654807 7.296408 6.886108 6.644738 7.423452

write.csv(rld_mat_heart_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/rlogtransformed_counts_PB.csv"), quote = FALSE, row.names = T)
rld_cor_heart_PB <- cor(rld_mat_heart_PB)
write.csv(rld_cor_heart_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/deseq2_pairwise_correlations_PB.csv"), quote = FALSE, row.names = T)
        
# Plot heatmap
png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_heatmap.png"), height=6, width=7.5, units="in", res=300)
pheatmap(rld_cor_heart_PB, annotation = cluster_metadata_heart_PB[, c("sex_id"), drop=F])
dev.off()
        
# Run DESeq2 differential expression analysis
dds_deseq_heart_PB <- DESeq(dds_heart_PB)
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
        
# Plot dispersion estimates
png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_dispersion_plot.png"))
plotDispEsts(dds_deseq_heart_PB)
dev.off()

# Output results of Wald test for contrast for A vs B
contrast_heart_PB <- c("sex_id", levels(cluster_metadata_heart_PB$sex_id)[2], levels(cluster_metadata_heart_PB$sex_id)[1])
        [1] "sex_id" "Male"   "Female"
# resultsNames(dds)
res_heart_PB <- results(dds_deseq_heart_PB, 
                contrast = contrast_heart_PB,
                alpha = 0.05)

res_heart_PB_ashr <- lfcShrink(dds_deseq_heart_PB, 
                contrast=contrast_heart_PB,
                res=res_heart_PB, type="ashr")
                     
 # Try other adjsuting methods
 res_heart_PB_bon <- results(dds_deseq_heart_PB, 
                contrast = contrast_heart_PB,
                alpha = 0.05,
                pAdjustMethod="bonferroni")                   
                     
                      
# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1
        
# Turn the results object into a tibble for use with tidyverse functions

res_tbl_heart_PB <- res_heart_PB %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
write.csv(res_tbl_heart_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_all_genes.csv"), quote = FALSE, row.names = FALSE)


res_tbl_heart_PB_ashr <- res_heart_PB_ashr %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
write.csv(res_tbl_heart_PB_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_all_genes.csv"), quote = FALSE, row.names = FALSE)
        
# Subset the significant results
sig_res_heart_PB <- dplyr::filter(res_tbl_heart_PB, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_heart_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes.csv"), quote = FALSE, row.names = FALSE)

sig_res_heart_PB_ashr <- dplyr::filter(res_tbl_heart_PB_ashr, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_heart_PB_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_heart <- dplyr::filter(sig_res_heart_PB, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_heart <- dplyr::filter(sig_res_heart_PB, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_ashr_heart <- dplyr::filter(sig_res_heart_PB_ashr, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_ashr_heart <- dplyr::filter(sig_res_heart_PB_ashr, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_ashr_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_ashr_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

## ggplot of top genes
normalized_counts_heart_PB <- counts(dds_deseq_heart_PB, normalized = TRUE)
write.csv(normalized_counts_heart_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/dds_normalized_unlogged_counts_PB.csv"), quote = FALSE, row.names = T)       
    
## Order results by padj values
top20_sig_genes_heart_PB_ashr <- sig_res_heart_PB_ashr %>%
                dplyr::arrange(padj) %>%
                dplyr::pull(gene) %>%
                head(n=20)

top20_sig_norm_heart_PB_ashr <- data.frame(normalized_counts_heart_PB) %>%
                rownames_to_column(var = "gene") %>%
                dplyr::filter(gene %in% top20_sig_genes_heart_PB_ashr)
        
gathered_top20_sig_heart_PB_ashr <- top20_sig_norm_heart_PB_ashr %>%
                gather(colnames(top20_sig_norm_heart_PB_ashr)[2:length(colnames(top20_sig_norm_heart_PB_ashr))], key = "samplename", value = "normalized_counts")
        
gathered_top20_sig_heart_PB_ashr <- inner_join(ei_heart_PB[, c("sample_id", "sex_id" )], gathered_top20_sig_heart_PB_ashr, by = c("sample_id" = "samplename"))
        
## plot using ggplot2    
ggplot(gathered_top20_sig_heart_PB_ashr) +
                geom_point(aes(x = gene, 
                     y = normalized_counts, 
                     color = sex_id), 
                     position=position_jitter(w=0.1,h=0)) +
                scale_y_log10() +
                xlab("Genes") +
                ylab("log10 Normalized Counts") +
                ggtitle("Top 20 Significant DE Genes") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_shrinkashr_top20_DE_genes.png"))  

        str("volcano plot")
		res_tbl_thresh <- res_tbl_heart_PB_ashr[!is.na(res_tbl_heart_PB_ashr$padj), ] %>% mutate(threshold=padj<0.05 & abs(log2FoldChange) >=1)
		min(log10(res_tbl_thresh$padj))
		
		ggplot(res_tbl_thresh) +
			geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
			ggtitle("Volcano plot") +
			xlab("log2 fold change") +
			ylab("-log10 adjusted p-value")+
			scale_color_manual(values=c("grey60", "red3")) +
			theme(legend.position="none", plot.title=element_text(size=rel(1.3), hjust=0.5), axis.title=element_text(size=rel(1.15)))
		ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_volcano_plot.png"))



