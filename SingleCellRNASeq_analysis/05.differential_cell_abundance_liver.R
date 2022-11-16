
liver_umap_rename_clusters_dbrmv_celltypes <- liver_umap_rename_clusters_dbrmv
liver_umap_rename_clusters_dbrmv_celltypes$CellType <- NA
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "0"))] <- "Hepatocyte 1"  
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "1"))] <- "Hepatocyte 1" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "2"))] <- "T lymphocyte" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "3"))] <- "Hepatocyte 2" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "4"))] <- "Hepatocyte 2" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "5"))] <- "Hepatocyte 2" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "6"))] <- "Macrophage" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "7"))] <- "B lymphocyte" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "8"))] <- "Biliary epithelial" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "9"))] <- "Erythrocyte" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "10"))] <- "Granulocyte" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "11"))] <- "Hepatocyte 1" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "12"))] <- "Erythrocyte" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "13"))] <- "Endothelial" 
liver_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(liver_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.3, "14"))] <- "Endothelial" 

table(liver_umap_rename_clusters_dbrmv_celltypes$CellType)
      B lymphocyte Biliary epithelial        Endothelial        Erythrocyte        Granulocyte       Hepatocyte 1       Hepatocyte 2         Macrophage       T lymphocyte 
              1509               1070                315               1191                636              17174               7939               1440               4496 
              
table(liver_umap_rename_clusters_dbrmv_celltypes$sample)
LiverF1 LiverF2 LiverF3 LiverM1 LiverM2 LiverM3 
   5814    5991    6075    4985    6070    6835 


## Proportions test (prop.test) or z-test method (similar to chisq.test) - Females:Males

tab_liver_Blymph <- matrix(c(682, 827, 17198, 17063), ncol=2)
data:  tab_liver_Blymph
X-squared = 14.263, df = 1, p-value = 0.000159
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.012304970 -0.003862561
sample estimates:
    prop 1     prop 2 
0.03814318 0.04622694 

tab_liver_bilep <- matrix(c(485, 585, 17395, 17305), ncol=2)
prop.test(tab_liver_bilep)
data:  tab_liver_bilep
X-squared = 9.3853, df = 1, p-value = 0.002187
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.009160567 -0.001988538
sample estimates:
    prop 1     prop 2 
0.02712528 0.03269983 

tab_liver_endothelial <- matrix(c(93, 222, 17787, 17668), ncol=2)
prop.test(tab_liver_endothelial)
data:  tab_liver_endothelial
X-squared = 52.403, df = 1, p-value = 4.521e-13
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.009198471 -0.005217179
sample estimates:
     prop 1      prop 2 
0.005201342 0.012409167 

tab_liver_eryt <- matrix(c(762, 429, 17118, 17461), ncol=2)
prop.test(tab_liver_eryt)
data:  tab_liver_eryt
X-squared = 95.927, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.01486795 0.02240720
sample estimates:
    prop 1     prop 2 
0.04261745 0.02397988 

tab_liver_hep1 <- matrix(c(8870, 8304, 9010, 9586), ncol=2)
prop.test(tab_liver_hep1)
data:  tab_liver_hep1
X-squared = 36.364, df = 1, p-value = 1.637e-09
alternative hypothesis: two.sided
95 percent confidence interval:
 0.02150956 0.04232060
sample estimates:
   prop 1    prop 2 
0.4960850 0.4641699 
		
tab_liver_hep2 <- matrix(c(3263, 4676, 14617, 13214), ncol=2)
prop.test(tab_liver_hep2)	
data:  tab_liver_hep2
X-squared = 321.76, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.08751024 -0.07025108
sample estimates:
   prop 1    prop 2 
0.1824944 0.2613751 

tab_liver_macro <- matrix(c(612, 828, 17268, 17062), ncol=2)
prop.test(tab_liver_macro)	
data:  tab_liver_macro
X-squared = 33.322, df = 1, p-value = 7.809e-09
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.016182455 -0.007926848
sample estimates:
    prop 1     prop 2 
0.03422819 0.04628284 

tab_liver_granu <- matrix(c(300,336, 17580, 17554), ncol=2)
prop.test(tab_liver_granu)		
data:  tab_liver_granu
X-squared = 1.9411, df = 1, p-value = 0.1635
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.0047977105  0.0007918732
sample estimates:
    prop 1     prop 2 
0.01677852 0.01878144 

tab_liver_tlymph <- matrix(c(2813, 1683, 15067, 16207), ncol=2)
prop.test(tab_liver_tlymph)	
data:  tab_liver_tlymph
X-squared = 324.99, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.05635597 0.07014747
sample estimates:
   prop 1    prop 2 
0.1573266 0.0940749 

