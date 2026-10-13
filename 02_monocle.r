# required format
rm(list = ls())
gc()

ORIGINAL_DIR <- "/data/nas1/ouyangkang_OD/project/02_project/"
output <- file.path(ORIGINAL_DIR, '13_pseudotime')

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)


# code
library(Seurat)
library(monocle)
library(ggplot2)
library(cowplot)
library(Matrix)
set.seed(1026)

count_matrix <- readMM("./02_monocle_matrix.mtx")
barcode_info <- read.table('./02_monocle_meta.tsv.gz',sep = '\t', row.names = 1,
                           col.names = c('index', 'group', 'celltype', 'cnv_leiden', 'cnv_score'))
feature_info <- read.table('./02_monocle_feature.tsv.gz', sep = '\t', col.names = c('gene'))

rownames(count_matrix) <- feature_info$gene
colnames(count_matrix) <- rownames(barcode_info)

count_matrix <- as(count_matrix, 'CsparseMatrix')


pd <- new('AnnotatedDataFrame', data=barcode_info)
fd <- new('AnnotatedDataFrame',
          data=data.frame(gene_short_name = feature_info$gene,
                          gene_id = feature_info$gene,
                          row.names = feature_info$gene,
                          stringsAsFactors = FALSE))

cds <- newCellDataSet(count_matrix,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, min_cells_detected = 15)
disp_tbl <- dispersionTable(cds)
ordering_genes <- disp_tbl$gene_id[disp_tbl$mean_expression >= 0.15]
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)


# plot
p1 <- plot_cell_trajectory(cds, color_by = "group") +
  ggtitle("group") +
  scale_color_manual(values = c(
    "Cancer" = "#FF3300",
    "Normal" = "#57C3F3"))
p1
ggsave('./fig3_group_monocle.pdf', width = 6, height = 6)

p1 <- plot_cell_trajectory(cds, color_by = "celltype") +
  ggtitle("celltype cluster") +
  scale_color_manual(values = c(
    "Epi_Neural_pro" = "#FB8072",
    "Epi_Neuron_like" = "#f781bf",
    "Epi_Neuroendocrine" = "#7BAFDE"))
p1
ggsave('./fig3_celltype_monocle.pdf', width = 6, height = 6)


pData(cds)$cnv_leiden <- factor(pData(cds)$cnv_leiden,
                                levels = unique(pData(cds)$cnv_leiden))

p1 <- plot_cell_trajectory(cds, color_by = "cnv_leiden") +
  ggtitle("CNV leiden cluster")+
  scale_color_manual(values = c(
    '1' = '#FB8072',
    '2' = '#7BAFDE',
    '6' = '#f781bf'
  ))
p1
ggsave('./fig3_cnv_leiden_monocle.pdf', width = 6, height = 6)

p1 <- plot_cell_trajectory(cds, color_by = "cnv_score") +
  ggtitle("CNV score cluster") +
  scale_color_gradient(low = "#57C3F3", high = "#E7298A")
p1
ggsave('./fig3_cnv_score_monocle.pdf', width = 6, height = 6)

p1 <- plot_cell_trajectory(cds, color_by = "State") +
  ggtitle("monocle state")
p1
ggsave('./fig3_state_monocle.pdf', width = 6, height = 6)

p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  ggtitle("pseudotime") +
  scale_color_gradient(low = "#57C3F3", high = "#E7298A")
p1
ggsave('./fig3_pseudotime_monocle.pdf', width = 6, height = 6)


# plot gene heatmap
diff_test <- differentialGeneTest(cds[ordering_genes, ],
                                  fullModelFormulaStr = '~State + Pseudotime',
                                  reducedModelFormulaStr = '~Pseudotime')
sig_genes <- row.names(subset(diff_test, qval < 0.01))[1:100]


library(RColorBrewer)
library(pheatmap)
my_palette <- colorRampPalette(c("#57C3F3", "white", "#E7298A"))(64)

pdf("./fig3_gene_heatmap_branch_1.pdf", width = 5, height = 10)
plot_genes_branched_heatmap(cds[sig_genes, ],
                            branch_point = 1,
                            num_clusters = 3,
                            hmcols = my_palette,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()

pdf("./fig3_gene_heatmap_branch_2.pdf", width = 5, height = 10)
plot_genes_branched_heatmap(cds[sig_genes, ],
                            branch_point = 2,
                            num_clusters = 3,
                            hmcols = my_palette,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()


# plot single gene
interested_genes <- c('MAP1B', "SYT1", "LSAMP", "NRXN1", "ENO1")
cds_subset <- cds[interested_genes, ]
p <- plot_genes_branched_pseudotime(cds_subset, branch_point = 1, color_by = 'celltype',
                                    branch_labels = c("Branch 1", "Branch 2")) +
  scale_color_manual(values = c(
      "Epi_Neural_pro" = "#FB8072",
      "Epi_Neuron_like" = "#f781bf",
      "Epi_Neuroendocrine" = "#7BAFDE"))
p
ggsave('./fig3_single_genes_branch_1_pseudotime.pdf', p, width = 6, height = 8)

p <- plot_genes_branched_pseudotime(cds_subset, branch_point = 2, color_by = 'celltype',
                                    branch_labels = c("Branch 1", "Branch 2")) +
  scale_color_manual(values = c(
    "Epi_Neural_pro" = "#FB8072",
    "Epi_Neuron_like" = "#f781bf",
    "Epi_Neuroendocrine" = "#7BAFDE"))
p
ggsave('./fig3_single_genes_branch_2_pseudotime.pdf', p, width = 6, height = 8)



# GO KEGG
ht <- plot_genes_branched_heatmap(cds[sig_genes, ],
                            branch_point = 2,
                            num_clusters = 3,
                            use_gene_short_name = T,
                            show_rownames = T,
                            return_heatmap = TRUE)
module_list <- split(ht$ph_res$tree_row$labels, ht$annotation_row$Cluster)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

gene_entrez <- bitr(module_list[[1]], fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Hs.eg.db) %>% pull(ENTREZID)

go_bp <- enrichGO(gene          = gene_entrez,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)

pdf('./fig3_branch_2_module_1_GO.pdf', width = 6, height = 6)
barplot(go_bp)
dev.off()

gene_entrez <- bitr(module_list[[2]], fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Hs.eg.db) %>% pull(ENTREZID)

go_bp <- enrichGO(gene          = gene_entrez,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)

pdf('./fig3_branch_2_module_2_GO.pdf', width = 6, height = 6)
barplot(go_bp)
dev.off()

gene_entrez <- bitr(module_list[[3]], fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Hs.eg.db) %>% pull(ENTREZID)

go_bp <- enrichGO(gene          = gene_entrez,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)

pdf('./fig3_branch_2_module_3_GO.pdf', width = 6, height = 6)
barplot(go_bp)
dev.off()
