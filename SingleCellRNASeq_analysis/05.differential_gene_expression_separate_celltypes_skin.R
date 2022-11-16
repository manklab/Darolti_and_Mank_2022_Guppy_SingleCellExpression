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

load("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/04.remove_doublets/skin_filtered/skin_umap_3_rename_clusters_dbrmv.RData")

head(skin_umap_3_rename_clusters_dbrmv)

counts_skin <- skin_umap_3_rename_clusters_dbrmv@assays$RNA@counts
metadata_skin <- skin_umap_3_rename_clusters_dbrmv@meta.data
metadata_skin$CellType <- factor(skin_umap_3_rename_clusters_dbrmv@active.ident)

sce_skin <- SingleCellExperiment(assays=list(counts=counts_skin), colData=metadata_skin)
class: SingleCellExperiment 
dim: 21975 12106 
metadata(0):
assays(1): counts
rownames(21975): ENSPREG00000006268 arl13b ... guca1d ENSPREG00000006408
rowData names(0):
colnames(12106): F1_Lib19_AAACCCACAAACCATC-1 F1_Lib19_AAACCCACAATGCAGG-1 ... M3_Lib33_TTTGTTGCAAATGGCG-1 M3_Lib33_TTTGTTGCATGCCATA-1
colData names(27): seq_folder nUMI ... seurat_clusters CellType
reducedDimNames(0):
altExpNames(0):

assays(sce_skin)
		List of length 1
		names(1): counts
dim(counts(sce_skin))
		[1] 21975 12106
counts(sce_skin)[1:2, 1:2]
		2 x 2 sparse Matrix of class "dgCMatrix"
		                   F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1
		ENSPREG00000006268                          .                          .
		arl13b                                      .                          .
dim(colData(sce_skin))
		[1] 12106    27


# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_skin)$sample_id <- as.factor(colData(sce_skin)$sample)
colData(sce_skin)$sex_id <- as.factor(colData(sce_skin)$sex)
colData(sce_skin)$cluster_id <- as.factor(colData(sce_skin)$CellType)

kids_skin <- purrr::set_names(levels(sce_skin$cluster_id))
kids_skin
             T lymphocyte                Fibroblast                   Stromal        Mitochondrial rich                 Epidermal              Keratinoycte       Mesenchymal stromal                Macrophage                Melanocyte 
           "T lymphocyte"              "Fibroblast"                 "Stromal"      "Mitochondrial rich"               "Epidermal"            "Keratinoycte"     "Mesenchymal stromal"              "Macrophage"              "Melanocyte" 
              Granulocyte              B lymphocyte   Epidermal (myelin-rich)           Skeletal muscle 
            "Granulocyte"            "B lymphocyte" "Epidermal (myelin-rich)"         "Skeletal muscle" 

nk_skin <- length(kids_skin)
nk_skin
		[1] 13

sids_skin <- purrr::set_names(levels(sce_skin$sample_id))
sids_skin
		   SkinF1   SkinF2   SkinF3   SkinM1   SkinM2   SkinM3 
		 "SkinF1" "SkinF2" "SkinF3" "SkinM1" "SkinM2" "SkinM3" 
ns_skin <- length(sids_skin)
ns_skin
		[1] 6

# Number of cells per sample
table(sce_skin$sample_id)
SkinF1 SkinF2 SkinF3 SkinM1 SkinM2 SkinM3 
  3193   1491   3747    917    503   2255 
  
# Turn named vector into a numeric vector
n_cells_skin <- as.numeric(table(sce_skin$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_skin <- match(sids_skin, sce_skin$sample_id)
m_skin
[1]    1 3194 4685 8432 9349 9852

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_skin <- data.frame(colData(sce_skin)[m_skin, ], n_cells_skin, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
dim(sce_skin)
		[1] 21975 12106
sce_skin <- sce_skin[rowSums(counts(sce_skin)>0)>0,]
dim(sce_skin)	
		[1] 18844 12106
sce_skin <- sce_skin[rowSums(counts(sce_skin)>1)>=10,]
dim(sce_skin)
		[1]  9297 12114
		

# Aggregate counts per sample_id and cluster_id
groups_skin <- colData(sce_skin)[, c("cluster_id", "sample_id")]
pb_skin <- aggregate.Matrix(t(counts(sce_skin)), groupings=groups_skin, fun="sum")
class(pb_skin)
		[1] "dgCMatrix"
		attr(,"package")
		[1] "Matrix"
dim(pb_skin)
		[1]    78 9297

pb_skin[1:6, 1:6]
		6 x 6 sparse Matrix of class "dgCMatrix"
                          zgc:152904 gart abcc4 timmdc1 VANGL1 sap18
T lymphocyte_SkinF1                .    1     1       4      .    20
Fibroblast_SkinF1                 12    2     2       6     72    33
Stromal_SkinF1                    10    3     1       5      9    27
Mitochondrial rich_SkinF1          3    2     1       .      2     6
Epidermal_SkinF1                   .    2     2       3      2    14
Keratinoycte_SkinF1                5    .     .       .      9     3


## Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_skin <- sapply(stringr::str_split(rownames(pb_skin), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_skin <- split.data.frame(pb_skin, factor(splitf_skin)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
class(pb_skin)
		[1] "list"
str(pb_skin)

# Number of cells for each cell type and sample
options(width=100)
table(sce_skin$cluster_id, sce_skin$sample_id)
                          SkinF1 SkinF2 SkinF3 SkinM1 SkinM2 SkinM3
  T lymphocyte               448    218   1147     72     91    363
  Fibroblast                 744    381    882    278    155    679
  Stromal                    358    172    374    169     57    244
  Mitochondrial rich         414    101    223     53     33    174
  Epidermal                  125    147    314     67     35    156
  Keratinoycte               347    183    113     40     42    105
  Mesenchymal stromal        407    105    176     54      6     63
  Macrophage                 110     73    152     48     16    116
  Melanocyte                  71     26     50     61     31    198
  Granulocyte                 43     38     99     21     13     48
  B lymphocyte                40     12     94     24      7     62
  Epidermal (myelin-rich)     72     24     37     23      4     18
  Skeletal muscle             14     11     86      7     13     29


get_sample_ids_skin <- function(x){
	pb_skin[[x]] %>% colnames()
}
de_samples_skin <- map(1:length(kids_skin), get_sample_ids_skin) %>% unlist()

samples_list_skin <- map(1:length(kids_skin), get_sample_ids_skin)

get_cluster_ids_skin <- function(x){
	rep(names(pb_skin)[x], each=length(samples_list_skin[[x]]))
}

de_cluster_ids_skin <- map(1:length(kids_skin), get_cluster_ids_skin) %>% unlist()

gg_df_skin <- data.frame(cluster_id=de_cluster_ids_skin, sample_id=de_samples_skin)
gg_df_skin <- left_join(gg_df_skin, ei_skin[, c("sample_id", "sex_id", "n_cells_skin")])
metadata_skin <- gg_df_skin %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_skin)

metadata_skin$cluster_id <- as.factor(metadata_skin$cluster_id)

clusters_skin <- levels(metadata_skin$cluster_id)

get_dds_resultsAvsB <- function(x, A, B){
        cluster_metadata_skin <- metadata_skin[which(metadata_skin$cluster_id == clusters_skin[x]), ]
        rownames(cluster_metadata_skin) <- cluster_metadata_skin$sample_id
        counts_skin <- pb_skin[[clusters_skin[x]]]
        cluster_counts_skin <- data.frame(counts_skin[, which(colnames(counts_skin) %in% rownames(cluster_metadata_skin))])
        
        all(rownames(cluster_metadata_skin) == colnames(cluster_counts_skin))        
        
        str("Create DESeq2 object")
        dds_skin <- DESeqDataSetFromMatrix(cluster_counts_skin, 
                                      colData = cluster_metadata_skin, 
                                      design = ~ sex_id)
        raw_counts <- counts(dds_skin)
        write.csv(raw_counts, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/",clusters_skin[x],"_dds_raw_counts.csv"), quote = FALSE, row.names = T)
        
        fpm_skin <- fpm(dds_skin)
        write.csv(fpm_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_dds_fpm.csv"), quote = FALSE, row.names = T)
        
        cpm_skin <- calculateCPM(dds_skin)
        write.csv(cpm_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_dds_cpm.csv"), quote = FALSE, row.names = T)
        
        
        str("Transform counts for data visualization")
        rld_skin <- rlog(dds_skin, blind=TRUE)
        
        str("Plot PCA")
        DESeq2::plotPCA(rld_skin, ntop=100, intgroup = "sex_id")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_skin[x], "_specific_PCAplot.png"))
        
        DESeq2::plotPCA(rld_skin, ntop=100, intgroup = "n_cells_skin")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_skin[x], "_specific_PCAplot_cellcount.png"))
        
        str("Extract the rlog matrix from the object and compute pairwise correlation values")
        rld_mat_skin <- assay(rld_skin)
        write.csv(rld_mat_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_rlogtransformed_counts.csv"), quote = FALSE, row.names = T)
        rld_cor_skin <- cor(rld_mat_skin)
        write.csv(rld_cor_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_pairwise_correlations.csv"), quote = FALSE, row.names = T)
        
        
        str("Plot heatmap")
        png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_skin[x], "_specific_heatmap_correl.png"), height=6, width=7.5, units="in", res=300)
        pheatmap(rld_cor_skin, annotation = cluster_metadata_skin[, c("sex_id"), drop=F])
        dev.off()

        str("Run DESeq2 differential expression analysis")
        dds_deseq_skin <- DESeq(dds_skin)
        
        str("Plot dispersion estimates")
        png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_skin[x], "_dispersion_plot.png"), height=4, width=6, units="in", res=300)
        plotDispEsts(dds_deseq_skin)
        dev.off()

        str("Output results of Wald test for contrast for A vs B")
        contrast_skin <- c("sex_id", levels(cluster_metadata_skin$sex_id)[A], levels(cluster_metadata_skin$sex_id)[B])
        str(contrast_skin)
        
        str("resultsNames(dds)")
        res_skin <- results(dds_deseq_skin, 
                       contrast = contrast_skin,
                       alpha = 0.05)
        
        res_skin_ashr <- lfcShrink(dds_deseq_skin, 
                         contrast=contrast_skin,
                         res=res_skin, type="ashr")
        str("Set thresholds")
        padj_cutoff <- 0.05
        log2fc_cutoff <- 1
        
        str("Turn the results object into a tibble for use with tidyverse functions")
        res_tbl_skin <- res_skin %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()
         write.csv(res_tbl_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_deseq_beforeshrink_all_genes.csv"), quote = FALSE, row.names = FALSE)
         
        res_tbl_skin_ashr <- res_skin_ashr %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()
        write.csv(res_tbl_skin_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_", levels(cluster_metadata_skin$sex_id)[2], "_vs_", levels(cluster_metadata_skin$sex_id)[1], "_shrinkashr_all_genes.csv"), quote = FALSE, row.names = FALSE)
        
        str("Subset the significant results")
        sig_res_skin <- dplyr::filter(res_tbl_skin, padj < padj_cutoff) %>% dplyr::arrange(padj) 
        write.csv(sig_res_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_deseq_beforeshrink_sig_genes.csv"), quote = FALSE, row.names = FALSE)
        
        sig_res_skin_ashr <- dplyr::filter(res_tbl_skin_ashr, padj < padj_cutoff) %>% dplyr::arrange(padj)
        write.csv(sig_res_skin_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_", levels(cluster_metadata_skin$sex_id)[A], "_vs_", levels(cluster_metadata_skin$sex_id)[B], "_shrinkashr_sig_genes.csv"), quote = FALSE, row.names = FALSE)
        
        sig_up_fc <- dplyr::filter(sig_res_skin, log2FoldChange >= log2fc_cutoff)
sig_down_fc <- dplyr::filter(sig_res_skin, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_deseq_beforeshrink_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_deseq_beforeshrink_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_ashr <- dplyr::filter(sig_res_skin_ashr, log2FoldChange >= log2fc_cutoff)
sig_down_fc_ashr <- dplyr::filter(sig_res_skin_ashr, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_deseq_shrinkashr_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_deseq_shrinkashr_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)
        
        str("ggplot of top genes")
        normalized_counts_skin <- counts(dds_deseq_skin, normalized = TRUE)
        write.csv(normalized_counts_skin, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_skin[x], "_rnormalized_counts.csv"), quote = FALSE, row.names = T)
        
        
        str("Order results by padj values")
        top20_sig_genes_skin_ashr <- sig_res_skin_ashr %>%
                dplyr::arrange(padj) %>%
                dplyr::pull(gene) %>%
                head(n=20)
        
        top20_sig_norm_skin_ashr <- data.frame(normalized_counts_skin) %>%
                rownames_to_column(var = "gene") %>%
                dplyr::filter(gene %in% top20_sig_genes_skin_ashr)
        
        gathered_top20_sig_skin_ashr <- top20_sig_norm_skin_ashr %>%
                gather(colnames(top20_sig_norm_skin_ashr)[2:length(colnames(top20_sig_norm_skin_ashr))], key = "samplename", value = "normalized_counts")
        
        gathered_top20_sig_skin_ashr <- inner_join(ei_skin[, c("sample_id", "sex_id" )], gathered_top20_sig_skin_ashr, by = c("sample_id" = "samplename"))
        
        str("plot using ggplot2")    
         ggplot(gathered_top20_sig_skin_ashr) +
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
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_skin[x], "_", levels(cluster_metadata_skin$sex_id)[A], "_vs_", levels(cluster_metadata_skin$sex_id)[B], "_shrinkashr_top20_DE_genes.png"))
		   
        str("volcano plot")
		res_tbl_thresh <- res_tbl_skin_ashr[!is.na(res_tbl_skin_ashr$padj), ] %>% mutate(threshold=padj<0.05 & abs(log2FoldChange) >= 1)
		min(log10(res_tbl_thresh$padj))
		
		ggplot(res_tbl_thresh) +
			geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
			ggtitle("Volcano plot") +
			xlab("log2 fold change") +
			ylab("-log10 adjusted p-value")+
			scale_color_manual(values=c("grey60", "red3")) +
			theme(legend.position="none", plot.title=element_text(size=rel(1.3), hjust=0.5), axis.title=element_text(size=rel(1.15)))
		ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/skin_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_skin[x], "_volcano_plot.png"))
		   

		
}      

map(1:length(clusters_skin), get_dds_resultsAvsB, A = 2, B = 1)






