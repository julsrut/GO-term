library(readr)
library(org.Hs.eg.db)
library(BiocManager)
library(AnnotationDbi)
BiocManager::install("clusterProfiler")
library(clusterProfiler)

Combo.df <- as.data.frame(Combo_truncated_DEseq2_Output)
gene <- names(Combo.df$ensLookup)
gene.df <- bitr(Combo.df, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
View(gene)
ego <- enrichGO(gene = gene, 
                universe =names(Combo.df),
                OrgDb = org.Hs.eg.db,
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)

#didnt work
ego <- enrichGO(gene          = gene,
                universe      = names(Combo.df),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego))
#worked!
ego2 <- enrichGO(gene         = Combo.df$ensembl_gene_id,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))
#also worked!
ego3 <- enrichGO(gene         = Combo.df$externalgenename,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego3))


#dope
dotplot(ego2, showCategory=30,
       font.size = 5)

#didn't work
enrichMap(ego2, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)

#didn't work
cnetplot(ego2, foldChange=Combo.df$`log2(FC)Combo`)

#didn't work
plotGOgraph(ego)

