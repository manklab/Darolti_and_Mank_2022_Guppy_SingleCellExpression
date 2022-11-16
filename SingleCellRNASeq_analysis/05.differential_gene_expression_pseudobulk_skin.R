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

load("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/04.remove_doublets/skin_filtered/skin_umap_3_rename_clusters_dbrmv.RData")

skin_umap_3_rename_clusters_dbrmv_PB <- skin_umap_3_rename_clusters_dbrmv
skin_umap_3_rename_clusters_dbrmv_PB$CellType <- "All"

head(skin_umap_3_rename_clusters_dbrmv_PB)

counts_skin_PB <- skin_umap_3_rename_clusters_dbrmv_PB@assays$RNA@counts
metadata_skin_PB <- skin_umap_3_rename_clusters_dbrmv_PB@meta.data
metadata_skin_PB$CellType <- "All"

sce_skin_PB <- SingleCellExperiment(assays=list(counts=counts_skin_PB), colData=metadata_skin_PB)
sce_skin_PB
class: SingleCellExperiment 
dim: 21975 12106 
metadata(0):
assays(1): counts
rownames(21975): ENSPREG00000006268 arl13b ... guca1d ENSPREG00000006408
rowData names(0):
colnames(12106): F1_Lib19_AAACCCACAAACCATC-1 F1_Lib19_AAACCCACAATGCAGG-1 ...
  M3_Lib33_TTTGTTGCAAATGGCG-1 M3_Lib33_TTTGTTGCATGCCATA-1
colData names(27): seq_folder nUMI ... seurat_clusters CellType
reducedDimNames(0):
altExpNames(0):

assays(sce_skin_PB)
		List of length 1
		names(1): counts
dim(counts(sce_skin_PB))
[1] 21975 12106
counts(sce_skin_PB)[1:2, 1:2]
		2 x 2 sparse Matrix of class "dgCMatrix"
		                   F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1
		ENSPREG00000006268                          .                          .
		arl13b                                      .                          .
dim(colData(sce_skin_PB))
[1] 12106    27

# Determine the number of clusters and the cluster names
colData(sce_skin_PB)$sample_id <- as.factor(colData(sce_skin_PB)$sample)
colData(sce_skin_PB)$sex_id <- as.factor(colData(sce_skin_PB)$sex)
colData(sce_skin_PB)$cluster_id <- as.factor(colData(sce_skin_PB)$CellType)
kids_skin_PB <- purrr::set_names(levels(sce_skin_PB$cluster_id))
kids_skin_PB
		  All 
		"All" 
nk_skin_PB <- length(kids_skin_PB)
nk_skin_PB
		[1] 1

sids_skin_PB <- purrr::set_names(levels(sce_skin_PB$sample_id))
sids_skin_PB
		  SkinF1   SkinF2   SkinF3   SkinM1   SkinM2   SkinM3 
		"SkinF1" "SkinF2" "SkinF3" "SkinM1" "SkinM2" "SkinM3" 
ns_skin_PB <- length(sids_skin_PB)
ns_skin_PB
		[1] 6

# Determine no. cells per sample
table(sce_skin_PB$sample_id)
		SkinF1 SkinF2 SkinF3 SkinM1 SkinM2 SkinM3 
  3193   1491   3747    917    503   2255

# Turn named vector into a numeric vector
n_cells_skin_PB <- as.numeric(table(sce_skin_PB$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_skin_PB <- match(sids_skin_PB, sce_skin_PB$sample_id)
m_skin_PB
[1]    1 3194 4685 8432 9349 9852

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_skin_PB <- data.frame(colData(sce_skin_PB)[m_skin_PB, ], n_cells_skin_PB, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
dim(sce_skin_PB)
[1] 21975 12106
sce_skin_PB <- sce_skin_PB[rowSums(counts(sce_skin_PB)>1)>=10,]
dim(sce_skin_PB)
[1]  9297 12106

# Aggregate counts per sample_id and cluster_id
groups_skin_PB <- colData(sce_skin_PB)[, c("cluster_id", "sample_id")]
pb_skin_PB <- aggregate.Matrix(t(counts(sce_skin_PB)), groupings=groups_skin_PB, fun="sum")
class(pb_skin_PB)
		[1] "dgCMatrix"
		attr(,"package")
		[1] "Matrix"
dim(pb_skin_PB)
[1]    6 9297

> pb_skin_PB[1:6, 1:6]
6 x 6 sparse Matrix of class "dgCMatrix"
           zgc:152904 gart abcc4 timmdc1 VANGL1 sap18
All_SkinF1         51   28    13      27     98   165
All_SkinF2         49   19    17      22     44   132
All_SkinF3         47   38    33      64     43   295
All_SkinM1         11   27    13      20     46    82
All_SkinM2          .   13     4      11     21    46
All_SkinM3         21   46    20      35     23   214


## Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_skin_PB <- sapply(stringr::str_split(rownames(pb_skin_PB), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_skin_PB <- split.data.frame(pb_skin_PB, factor(splitf_skin_PB)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb_skin_PB)
		[1] "list"
str(pb_skin_PB)

options(width=100)
table(sce_skin_PB$cluster_id, sce_skin_PB$sample_id)
		      SkinF1 SkinF2 SkinF3 SkinM1 SkinM2 SkinM3
  All   3193   1491   3747    917    503   2255
  
get_sample_ids_skin_PB <- function(x){
	pb_skin_PB[[x]] %>% colnames()
}
de_samples_skin_PB <- map(1:length(kids_skin_PB), get_sample_ids_skin_PB) %>% unlist()

samples_list_skin_PB <- map(1:length(kids_skin_PB), get_sample_ids_skin_PB)

get_cluster_ids_skin_PB <- function(x){
	rep(names(pb_skin_PB)[x], each=length(samples_list_skin_PB[[x]]))
}

de_cluster_ids_skin_PB <- map(1:length(kids_skin_PB), get_cluster_ids_skin_PB) %>% unlist()

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
write.csv(raw_counts_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/dds_raw_counts_PB.csv"), quote = FALSE, row.names = T)

#-----------
fpm_PB <- fpm(dds_skin_PB)
head(fpm_PB)
              SkinF1    SkinF2    SkinF3    SkinM1    SkinM2    SkinM3
zgc:152904  8.795731 14.622938  6.205859  4.993184  0.000000  4.265467
gart        4.829029  5.670119  5.017503 12.255998  9.469009  9.343403
abcc4       2.242049  5.073264  4.357305  5.901036  2.913541  4.062349
timmdc1     4.656563  6.565401  8.450531  9.078517  8.012239  7.109111
VANGL1     16.901600 13.130801  5.677700 20.880589 15.296092  4.671702
sap18      28.456776 39.392403 38.951665 37.221920 33.505725 43.467137

cpm_PB <- calculateCPM(dds_skin_PB)
head(cpm_PB)
              SkinF1    SkinF2    SkinF3    SkinM1    SkinM2    SkinM3
zgc:152904  8.071300 16.291762  5.351166  5.941504  0.000000  4.392294
gart        4.431302  6.317214  4.326475 14.583692  9.259266  9.621215
abcc4       2.057390  5.652244  3.757202  7.021778  2.849005  4.183137
timmdc1     4.273041  7.314669  7.286695 10.802735  7.834763  7.320490
VANGL1     15.509556 14.629337  4.895748 24.846291 14.957276  4.810607
sap18      26.113028 43.888012 33.587108 44.291214 32.763556 44.759565

write.csv(cpm_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/dds_cpm_PB.csv"), quote = FALSE, row.names = T)
#--------------


# Transform counts for data visualization
rld_skin_PB <- rlog(dds_skin_PB, blind=TRUE)
colData(rld_skin_PB)
DataFrame with 6 rows and 5 columns
       cluster_id   sample_id   sex_id n_cells_skin_PB sizeFactor
         <factor> <character> <factor>       <numeric>  <numeric>
SkinF1        All      SkinF1   Female            3193   1.623613
SkinF2        All      SkinF2   Female            1491   0.938308
SkinF3        All      SkinF3   Female            3747   2.120704
SkinM1        All      SkinM1   Male               917   0.616878
SkinM2        All      SkinM2   Male               503   0.384435
SkinM3        All      SkinM3   Male              2255   1.378596
        
# Plot PCA
DESeq2::plotPCA(rld_skin_PB, ntop=100,intgroup = "sex_id")
ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_PCAplot.png"))

        
        DESeq2::plotPCA(rld_skin_PB, ntop=100, intgroup = "n_cells_skin_PB")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_PCAplot_cellcount.png"))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat_skin_PB <- assay(rld_skin_PB)
head(rld_mat_skin_PB)
            SkinF1   SkinF2   SkinF3   SkinM1   SkinM2   SkinM3
zgc:152904 4.631737 4.962156 4.429290 4.324958 3.835573 4.244503
gart       4.500724 4.585374 4.516004 5.023624 4.853964 4.862439
abcc4      3.630193 3.942278 3.876991 4.006756 3.740410 3.845762
timmdc1    4.461708 4.634565 4.777100 4.810497 4.736546 4.676384
VANGL1     5.640770 5.455581 4.942916 5.787201 5.549363 4.859474
sap18      6.827956 7.090540 7.081900 7.043657 6.961250 7.173969

write.csv(rld_mat_skin_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/rlogtransformed_counts_PB.csv"), quote = FALSE, row.names = T)
rld_cor_skin_PB <- cor(rld_mat_skin_PB)
write.csv(rld_cor_skin_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/deseq2_pairwise_correlations_PB.csv"), quote = FALSE, row.names = T)
        
# Plot heatmap
png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_heatmap.png"), height=6, width=7.5, units="in", res=300)
pheatmap(rld_cor_skin_PB, annotation = cluster_metadata_skin_PB[, c("sex_id"), drop=F])
dev.off()
        
# Run DESeq2 differential expression analysis
dds_deseq_skin_PB <- DESeq(dds_skin_PB)
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
        
# Plot dispersion estimates
png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_dispersion_plot.png"))
plotDispEsts(dds_deseq_skin_PB)
dev.off()

# Output results of Wald test for contrast for A vs B
contrast_skin_PB <- c("sex_id", levels(cluster_metadata_skin_PB$sex_id)[2], levels(cluster_metadata_skin_PB$sex_id)[1])
        [1] "sex_id" "Male"   "Female"
# resultsNames(dds)
res_skin_PB <- results(dds_deseq_skin_PB, 
                contrast = contrast_skin_PB,
                alpha = 0.05)

res_skin_PB_ashr <- lfcShrink(dds_deseq_skin_PB, 
                contrast=contrast_skin_PB,
                res=res_skin_PB, type="ashr")
                                  
                                      
# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1
        
# Turn the results object into a tibble for use with tidyverse functions

res_tbl_skin_PB <- res_skin_PB %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
write.csv(res_tbl_skin_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_all_genes.csv"), quote = FALSE, row.names = FALSE)


res_tbl_skin_PB_ashr <- res_skin_PB_ashr %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
write.csv(res_tbl_skin_PB_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_all_genes.csv"), quote = FALSE, row.names = FALSE)
        
# Subset the significant results
sig_res_skin_PB <- dplyr::filter(res_tbl_skin_PB, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_skin_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes.csv"), quote = FALSE, row.names = FALSE)

sig_res_skin_PB_ashr <- dplyr::filter(res_tbl_skin_PB_ashr, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_skin_PB_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_skin <- dplyr::filter(sig_res_skin_PB, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_skin <- dplyr::filter(sig_res_skin_PB, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_ashr_skin <- dplyr::filter(sig_res_skin_PB_ashr, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_ashr_skin <- dplyr::filter(sig_res_skin_PB_ashr, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_ashr_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_ashr_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

## ggplot of top genes
normalized_counts_skin_PB <- counts(dds_deseq_skin_PB, normalized = TRUE)
write.csv(normalized_counts_skin_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/dds_normalized_unlogged_counts_PB.csv"), quote = FALSE, row.names = T)       
    
## Order results by padj values
top20_sig_genes_skin_PB_ashr <- sig_res_skin_PB_ashr %>%
                dplyr::arrange(padj) %>%
                dplyr::pull(gene) %>%
                head(n=20)

top20_sig_norm_skin_PB_ashr <- data.frame(normalized_counts_skin_PB) %>%
                rownames_to_column(var = "gene") %>%
                dplyr::filter(gene %in% top20_sig_genes_skin_PB_ashr)
        
gathered_top20_sig_skin_PB_ashr <- top20_sig_norm_skin_PB_ashr %>%
                gather(colnames(top20_sig_norm_skin_PB_ashr)[2:length(colnames(top20_sig_norm_skin_PB_ashr))], key = "samplename", value = "normalized_counts")
        
gathered_top20_sig_skin_PB_ashr <- inner_join(ei_skin_PB[, c("sample_id", "sex_id" )], gathered_top20_sig_skin_PB_ashr, by = c("sample_id" = "samplename"))
        
## plot using ggplot2    
ggplot(gathered_top20_sig_skin_PB_ashr) +
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
ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_shrinkashr_top20_DE_genes.png"))  

        str("volcano plot")
		res_tbl_thresh <- res_tbl_skin_PB_ashr[!is.na(res_tbl_skin_PB_ashr$padj), ] %>% mutate(threshold=padj<0.05 & abs(log2FoldChange) >=1)
		min(log10(res_tbl_thresh$padj))
		
		ggplot(res_tbl_thresh) +
			geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
			ggtitle("Volcano plot") +
			xlab("log2 fold change") +
			ylab("-log10 adjusted p-value")+
			scale_color_manual(values=c("grey60", "red3")) +
			theme(legend.position="none", plot.title=element_text(size=rel(1.3), hjust=0.5), axis.title=element_text(size=rel(1.15)))
		ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_volcano_plot.png"))



