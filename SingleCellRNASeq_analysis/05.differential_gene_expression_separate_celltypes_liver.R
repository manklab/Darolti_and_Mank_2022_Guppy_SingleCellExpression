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

load("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/04.remove_doublets/liver_filtered/liver_umap_3_rename_clusters_dbrmv.RData")

head(liver_umap_3_rename_clusters_dbrmv)

counts_liver <- liver_umap_3_rename_clusters_dbrmv@assays$RNA@counts
metadata_liver <- liver_umap_3_rename_clusters_dbrmv@meta.data
metadata_liver$CellType <- factor(liver_umap_3_rename_clusters_dbrmv@active.ident)

sce_liver <- SingleCellExperiment(assays=list(counts=counts_liver), colData=metadata_liver)
sce_liver 
class: SingleCellExperiment 
dim: 21975 35770 
metadata(0):
assays(1): counts
rownames(21975): ENSPREG00000006268 arl13b ... guca1d ENSPREG00000006408
rowData names(0):
colnames(35770): F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1 ...
  M3_Lib14_TTTGTTGTCGAGTCCG-1 M3_Lib14_TTTGTTGTCTTACACT-1
colData names(27): seq_folder nUMI ... seurat_clusters CellType
reducedDimNames(0):
altExpNames(0):

assays(sce_liver)
		List of length 1
		names(1): counts
dim(counts(sce_liver))
		[1] 21975 12106
counts(sce_liver)[1:2, 1:2]
		2 x 2 sparse Matrix of class "dgCMatrix"
		                   F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1
		ENSPREG00000006268                          .                          .
		arl13b                                      .                          .
dim(colData(sce_liver))
		[1] 12106    27


# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_liver)$sample_id <- as.factor(colData(sce_liver)$sample)
colData(sce_liver)$sex_id <- as.factor(colData(sce_liver)$sex)
colData(sce_liver)$cluster_id <- as.factor(colData(sce_liver)$CellType)

kids_liver <- purrr::set_names(levels(sce_liver$cluster_id))
kids_liver
        Hepatocyte 1         T lymphocyte         Hepatocyte 2           Macrophage 
      "Hepatocyte 1"       "T lymphocyte"       "Hepatocyte 2"         "Macrophage" 
        B lymphocyte   Biliary epithelial          Erythrocyte          Granulocyte 
      "B lymphocyte" "Biliary epithelial"        "Erythrocyte"        "Granulocyte" 
         Endothelial 
       "Endothelial" 

nk_liver <- length(kids_liver)
nk_liver
		[1] 9

sids_liver <- purrr::set_names(levels(sce_liver$sample_id))
sids_liver
		   LiverF1   LiverF2   LiverF3   LiverM1   LiverM2   LiverM3 
		 "LiverF1" "LiverF2" "LiverF3" "LiverM1" "LiverM2" "LiverM3" 
ns_liver <- length(sids_liver)
ns_liver
		[1] 6

# Number of cells per sample
table(sce_liver$sample_id)
LiverF1 LiverF2 LiverF3 LiverM1 LiverM2 LiverM3 
   5814    5991    6075    4985    6070    6835
  
# Turn named vector into a numeric vector
n_cells_liver <- as.numeric(table(sce_liver$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_liver <- match(sids_liver, sce_liver$sample_id)
m_liver
[1]     1  5815 11806 17881 22866 28936

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_liver <- data.frame(colData(sce_liver)[m_liver, ], n_cells_liver, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
dim(sce_liver)
[1] 21975 35770
sce_liver <- sce_liver[rowSums(counts(sce_liver)>0)>0,]
dim(sce_liver)	
[1] 18996 35770
sce_liver <- sce_liver[rowSums(counts(sce_liver)>1)>=10,]
dim(sce_liver)
[1] 11478 35770
		

# Aggregate counts per sample_id and cluster_id
groups_liver <- colData(sce_liver)[, c("cluster_id", "sample_id")]
pb_liver <- aggregate.Matrix(t(counts(sce_liver)), groupings=groups_liver, fun="sum")
class(pb_liver)
		[1] "dgCMatrix"
		attr(,"package")
		[1] "Matrix"
dim(pb_liver)
[1]    54 11478

pb_liver[1:6, 1:6]
		6 x 6 sparse Matrix of class "dgCMatrix"
                           ENSPREG00000006268 gart n6amt1 timmdc1 VANGL1 ENSPREG00000006674
Hepatocyte 1_LiverF1                        4    .      2      18      1                446
T lymphocyte_LiverF1                        .    3      7      12      4                135
Hepatocyte 2_LiverF1                        .    6      3      15      3               5464
Macrophage_LiverF1                          .    1      .      14      .                 43
B lymphocyte_LiverF1                        .    .      1       2      .                 20
Biliary epithelial_LiverF1                  .    1      .       5     17                 79


## Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_liver <- sapply(stringr::str_split(rownames(pb_liver), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_liver <- split.data.frame(pb_liver, factor(splitf_liver)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
class(pb_liver)
		[1] "list"
str(pb_liver)

# Number of cells for each cell type and sample
options(width=100)
table(sce_liver$cluster_id, sce_liver$sample_id)
                     LiverF1 LiverF2 LiverF3 LiverM1 LiverM2 LiverM3
  Hepatocyte 1          2963    3189    2718    3000    2364    2940
  T lymphocyte           997     778    1038     438     257     988
  Hepatocyte 2           822    1062    1379    1011    2407    1258
  Macrophage             223     204     185     199     194     435
  B lymphocyte           219     207     256      62     220     545
  Biliary epithelial     128     183     174     157     222     206
  Erythrocyte            322     233     207      21     279     129
  Granulocyte            114      99      87      50      85     201
  Endothelial             26      36      31      47      42     133


get_sample_ids_liver <- function(x){
	pb_liver[[x]] %>% colnames()
}
de_samples_liver <- map(1:length(kids_liver), get_sample_ids_liver) %>% unlist()

samples_list_liver <- map(1:length(kids_liver), get_sample_ids_liver)

get_cluster_ids_liver <- function(x){
	rep(names(pb_liver)[x], each=length(samples_list_liver[[x]]))
}

de_cluster_ids_liver <- map(1:length(kids_liver), get_cluster_ids_liver) %>% unlist()

gg_df_liver <- data.frame(cluster_id=de_cluster_ids_liver, sample_id=de_samples_liver)
gg_df_liver <- left_join(gg_df_liver, ei_liver[, c("sample_id", "sex_id", "n_cells_liver")])
metadata_liver <- gg_df_liver %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_liver)

metadata_liver$cluster_id <- as.factor(metadata_liver$cluster_id)

clusters_liver <- levels(metadata_liver$cluster_id)

get_dds_resultsAvsB <- function(x, A, B){
        cluster_metadata_liver <- metadata_liver[which(metadata_liver$cluster_id == clusters_liver[x]), ]
        rownames(cluster_metadata_liver) <- cluster_metadata_liver$sample_id
        counts_liver <- pb_liver[[clusters_liver[x]]]
        cluster_counts_liver <- data.frame(counts_liver[, which(colnames(counts_liver) %in% rownames(cluster_metadata_liver))])
        
        all(rownames(cluster_metadata_liver) == colnames(cluster_counts_liver))        
        
        str("Create DESeq2 object")
        dds_liver <- DESeqDataSetFromMatrix(cluster_counts_liver, 
                                      colData = cluster_metadata_liver, 
                                      design = ~ sex_id)
        raw_counts <- counts(dds_liver)
        write.csv(raw_counts, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/",clusters_liver[x],"_dds_raw_counts.csv"), quote = FALSE, row.names = T)
        
        fpm_liver <- fpm(dds_liver)
        write.csv(fpm_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_dds_fpm.csv"), quote = FALSE, row.names = T)
        
        cpm_liver <- calculateCPM(dds_liver)
        write.csv(cpm_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_dds_cpm.csv"), quote = FALSE, row.names = T)
        
        
        str("Transform counts for data visualization")
        rld_liver <- rlog(dds_liver, blind=TRUE)
        
        str("Plot PCA")
        DESeq2::plotPCA(rld_liver, ntop=100, intgroup = "sex_id")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_liver[x], "_specific_PCAplot.png"))
        
        DESeq2::plotPCA(rld_liver, ntop=100, intgroup = "n_cells_liver")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_liver[x], "_specific_PCAplot_cellcount.png"))
        
        str("Extract the rlog matrix from the object and compute pairwise correlation values")
        rld_mat_liver <- assay(rld_liver)
        write.csv(rld_mat_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_rlogtransformed_counts.csv"), quote = FALSE, row.names = T)
        rld_cor_liver <- cor(rld_mat_liver)
        write.csv(rld_cor_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_pairwise_correlations.csv"), quote = FALSE, row.names = T)
        
        
        str("Plot heatmap")
        png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_liver[x], "_specific_heatmap_correl.png"), height=6, width=7.5, units="in", res=300)
        pheatmap(rld_cor_liver, annotation = cluster_metadata_liver[, c("sex_id"), drop=F])
        dev.off()

        str("Run DESeq2 differential expression analysis")
        dds_deseq_liver <- DESeq(dds_liver)
        
        str("Plot dispersion estimates")
        png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_liver[x], "_dispersion_plot.png"), height=4, width=6, units="in", res=300)
        plotDispEsts(dds_deseq_liver)
        dev.off()

        str("Output results of Wald test for contrast for A vs B")
        contrast_liver <- c("sex_id", levels(cluster_metadata_liver$sex_id)[A], levels(cluster_metadata_liver$sex_id)[B])
        str(contrast_liver)
        
        str("resultsNames(dds)")
        res_liver <- results(dds_deseq_liver, 
                       contrast = contrast_liver,
                       alpha = 0.05)
        
        res_liver_ashr <- lfcShrink(dds_deseq_liver, 
                         contrast=contrast_liver,
                         res=res_liver, type="ashr")
        str("Set thresholds")
        padj_cutoff <- 0.05
        log2fc_cutoff <- 1
        
        str("Turn the results object into a tibble for use with tidyverse functions")
        res_tbl_liver <- res_liver %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()
         write.csv(res_tbl_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_deseq_beforeshrink_all_genes.csv"), quote = FALSE, row.names = FALSE)
         
        res_tbl_liver_ashr <- res_liver_ashr %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()
        write.csv(res_tbl_liver_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_", levels(cluster_metadata_liver$sex_id)[2], "_vs_", levels(cluster_metadata_liver$sex_id)[1], "_shrinkashr_all_genes.csv"), quote = FALSE, row.names = FALSE)
        
        str("Subset the significant results")
        sig_res_liver <- dplyr::filter(res_tbl_liver, padj < padj_cutoff) %>% dplyr::arrange(padj) 
        write.csv(sig_res_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_deseq_beforeshrink_sig_genes.csv"), quote = FALSE, row.names = FALSE)
        
        sig_res_liver_ashr <- dplyr::filter(res_tbl_liver_ashr, padj < padj_cutoff) %>% dplyr::arrange(padj)
        write.csv(sig_res_liver_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_", levels(cluster_metadata_liver$sex_id)[A], "_vs_", levels(cluster_metadata_liver$sex_id)[B], "_shrinkashr_sig_genes.csv"), quote = FALSE, row.names = FALSE)
        
        sig_up_fc <- dplyr::filter(sig_res_liver, log2FoldChange >= log2fc_cutoff)
sig_down_fc <- dplyr::filter(sig_res_liver, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_deseq_beforeshrink_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_deseq_beforeshrink_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_ashr <- dplyr::filter(sig_res_liver_ashr, log2FoldChange >= log2fc_cutoff)
sig_down_fc_ashr <- dplyr::filter(sig_res_liver_ashr, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_deseq_shrinkashr_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_deseq_shrinkashr_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)
        
        str("ggplot of top genes")
        normalized_counts_liver <- counts(dds_deseq_liver, normalized = TRUE)
        write.csv(normalized_counts_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_liver[x], "_rnormalized_counts.csv"), quote = FALSE, row.names = T)
        
        
        str("Order results by padj values")
        top20_sig_genes_liver_ashr <- sig_res_liver_ashr %>%
                dplyr::arrange(padj) %>%
                dplyr::pull(gene) %>%
                head(n=20)
        
        top20_sig_norm_liver_ashr <- data.frame(normalized_counts_liver) %>%
                rownames_to_column(var = "gene") %>%
                dplyr::filter(gene %in% top20_sig_genes_liver_ashr)
        
        gathered_top20_sig_liver_ashr <- top20_sig_norm_liver_ashr %>%
                gather(colnames(top20_sig_norm_liver_ashr)[2:length(colnames(top20_sig_norm_liver_ashr))], key = "samplename", value = "normalized_counts")
        
        gathered_top20_sig_liver_ashr <- inner_join(ei_liver[, c("sample_id", "sex_id" )], gathered_top20_sig_liver_ashr, by = c("sample_id" = "samplename"))
        
        str("plot using ggplot2")    
         ggplot(gathered_top20_sig_liver_ashr) +
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
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_liver[x], "_", levels(cluster_metadata_liver$sex_id)[A], "_vs_", levels(cluster_metadata_liver$sex_id)[B], "_shrinkashr_top20_DE_genes.png"))
		   
        str("volcano plot")
		res_tbl_thresh <- res_tbl_liver_ashr[!is.na(res_tbl_liver_ashr$padj), ] %>% mutate(threshold=padj<0.05 & abs(log2FoldChange) >= 1)
		min(log10(res_tbl_thresh$padj))
		
		ggplot(res_tbl_thresh) +
			geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
			ggtitle("Volcano plot") +
			xlab("log2 fold change") +
			ylab("-log10 adjusted p-value")+
			scale_color_manual(values=c("grey60", "red3")) +
			theme(legend.position="none", plot.title=element_text(size=rel(1.3), hjust=0.5), axis.title=element_text(size=rel(1.15)))
		ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_liver[x], "_volcano_plot.png"))
		   

		
}      

map(1:length(clusters_liver), get_dds_resultsAvsB, A = 2, B = 1)






