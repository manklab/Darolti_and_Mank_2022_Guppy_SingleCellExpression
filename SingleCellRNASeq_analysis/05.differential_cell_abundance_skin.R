
skin_umap_rename_clusters_dbrmv_celltypes <- skin_umap_rename_clusters_dbrmv
skin_umap_rename_clusters_dbrmv_celltypes$CellType <- NA
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "0"))] <- "T lymphocyte"  
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "1"))] <- "Fibroblast" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "2"))] <- "Stromal" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "3"))] <- "Fibroblast" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "4"))] <- "Mitochondrial-rich" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "5"))] <- "Epidermal" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "6"))] <- "Keratinocyte" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "7"))] <- "Mesenchymal stromal" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "8"))] <- "Macrophage" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "9"))] <- "Melanocyte" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "10"))] <- "T lymphocyte" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "11"))] <- "Granulocyte" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "12"))] <- "B lymphocyte" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "13"))] <- "Epidermal (myelin-rich)" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "14"))] <- "Skeletal muscle" 
skin_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(skin_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.4, "15"))] <- "Mesenchymal stromal" 

table(skin_umap_rename_clusters_dbrmv_celltypes$CellType)
           B lymphocyte               Epidermal Epidermal (myelin-rich)              Fibroblast             Granulocyte            Keratinocyte              Macrophage              Melanocyte 
                    239                     844                     178                    3119                     262                     830                     515                     437 
    Mesenchymal stromal      Mitochondrial-rich         Skeletal muscle                 Stromal            T lymphocyte 
                    811                     998                     160                    1374                    2339

table(skin_umap_rename_clusters_dbrmv_celltypes$sample)
SkinF1 SkinF2 SkinF3 SkinM1 SkinM2 SkinM3 
  3193   1491   3747    917    503   2255 


## Proportions test (prop.test) or z-test method (similar to chisq.test) - Females:Males

tab_skin_epidermal <- matrix(c(586, 258, 7845, 3417), ncol=2)
prop.test(tab_skin_epidermal)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_epidermal
		X-squared = 0.0099972, df = 1, p-value = 0.9204
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 -0.010778377  0.009381007
		sample estimates:
		    prop 1     prop 2 
		0.06950540 0.07020408

tab_skin_fibroblast <- matrix(c(2007, 1112, 6424, 2563), ncol=2)
prop.test(tab_skin_fibroblast)
data:  tab_skin_fibroblast
X-squared = 55.394, df = 1, p-value = 9.863e-14
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.08214384 -0.04692612
sample estimates:
   prop 1    prop 2 
0.2380501 0.3025850 

tab_skin_granulocyte <- matrix(c(180, 82, 8251, 3593), ncol=2)
prop.test(tab_skin_granulocyte)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_granulocyte
		X-squared = 0.071255, df = 1, p-value = 0.7895
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 -0.006843853  0.004917564
		sample estimates:
		    prop 1     prop 2 
		0.02134978 0.02231293

tab_skin_blymph <- matrix(c(146, 93, 8285, 3582), ncol=2)
prop.test(tab_skin_blymph)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_blymph
		X-squared = 8.0332, df = 1, p-value = 0.004593
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 -0.013975514 -0.002002642
		sample estimates:
		    prop 1     prop 2 
		0.01731704 0.02530612

tab_skin_keratinocyte <- matrix(c(643, 187, 7788, 3488), ncol=2)
prop.test(tab_skin_keratinocyte)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_keratinocyte
		X-squared = 25.424, df = 1, p-value = 4.602e-07
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 0.01609898 0.03466463
		sample estimates:
		    prop 1     prop 2 
		0.07626616 0.05088435

tab_skin_macrophage <- matrix(c(335, 180, 8096, 3495), ncol=2)
prop.test(tab_skin_macrophage)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_macrophage
		X-squared = 5.1463, df = 1, p-value = 0.0233
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 -0.0175693172 -0.0009212388
		sample estimates:
		    prop 1     prop 2 
		0.03973431 0.04897959 

tab_skin_melanocyte <- matrix(c(147, 290, 8284, 3385), ncol=2)
prop.test(tab_skin_melanocyte)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_melanocyte
		X-squared = 276.23, df = 1, p-value < 2.2e-16
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 -0.07082455 -0.05212727
		sample estimates:
		    prop 1     prop 2 
		0.01743565 0.07891156 
		
tab_skin_messtromal <- matrix(c(688, 123, 7743, 3552), ncol=2)
prop.test(tab_skin_messtromal)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_messtromal
		X-squared = 94.103, df = 1, p-value < 2.2e-16
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 0.03969497 0.05657347
		sample estimates:
		    prop 1     prop 2 
		0.08160361 0.03346939

tab_skin_mito <- matrix(c(738, 260, 7693, 3415), ncol=2)
prop.test(tab_skin_mito)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_mito
		X-squared = 9.3129, df = 1, p-value = 0.002275
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 0.006337968 0.027233635
		sample estimates:
		   prop 1    prop 2 
		0.0875341 0.0707483 

tab_skin_epidmyel <- matrix(c(133, 45, 8298, 3630), ncol=2)
prop.test(tab_skin_epidmyel)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_pigment
		X-squared = 1.9647, df = 1, p-value = 0.161
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 -0.001105537  0.008165972
		sample estimates:
		    prop 1     prop 2 
		0.01577512 0.01224490

tab_skin_ribo <- matrix(c(1813, 526, 6618, 3149), ncol=2)
prop.test(tab_skin_ribo)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_ribo
		X-squared = 84.444, df = 1, p-value < 2.2e-16
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 0.05739351 0.08642745
		sample estimates:
		   prop 1    prop 2 
		0.2150397 0.1431293

tab_skin_skeletal <- matrix(c(111, 49, 8320, 3626), ncol=2)
prop.test(tab_skin_skeletal)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_skeletal
		X-squared = 1.9891e-29, df = 1, p-value = 1
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 -0.004770494  0.004435224
		sample estimates:
		    prop 1     prop 2 
		0.01316570 0.01333333

tab_skin_stromal <- matrix(c(904, 470, 7527, 3205), ncol=2)
prop.test(tab_skin_stromal)
			2-sample test for equality of proportions with continuity correction
		data:  tab_skin_stromal
		X-squared = 10.661, df = 1, p-value = 0.001094
		alternative hypothesis: two.sided
		95 percent confidence interval:
		 -0.03352032 -0.00781531
		sample estimates:
		   prop 1    prop 2 
		0.1072233 0.1278912

