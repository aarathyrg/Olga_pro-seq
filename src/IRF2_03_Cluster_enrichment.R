library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(enrichplot)
library(here)
library(readr)
library(clusterProfiler)
library(enrichR)
library(circlize)
library(grid)
library(biomaRt)
source("src/Ag_optimized_theme.R")
InDir1 <- here("Results/IRF2_01_DE/Processed_data/")
InDir2 <- here("Results/IRF2_02_Heatmaps/")

outdir <- here("Results/IRF2_03_Cluster_enrichment/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

#Input datasets
dataVoom <- read_rds(paste0(InDir1,"/dataVoom.rds"))
genes_by_cluster <- read_rds(paste0(InDir2,"/genes_by_cluster.rds"))
ht_obj <- read_rds(paste0(InDir2,"/ht_object.rds"))


#relevant functions
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
convert_to_symbol <- function(genes) {
  
  for (i in 1:3) {
    res <- tryCatch({
      getBM(
        attributes = c("ensembl_gene_id", "mgi_symbol"),
        filters = "ensembl_gene_id",
        values = genes,
        mart = mart
      )
    }, error = function(e) NULL)
    
    if (!is.null(res)) {
      return(unique(res$mgi_symbol[res$mgi_symbol != ""]))
    }
    
    Sys.sleep(2)
  }
  
  return(character(0))
}
#enrichment cluster profiler-------------------
cluster_dir <- paste0(outdir,"/Heatmap_clustering", "6")
if (!dir.exists(cluster_dir)) dir.create(cluster_dir)
Enrichment_dir <- file.path(cluster_dir, "enrichment_cluster_profiler")
if (!dir.exists(Enrichment_dir)) dir.create(Enrichment_dir)
# All genes in your dataset
bg_genes <- rownames(dataVoom$E)  # Ensembl IDs

enrichment_methods <- list(
  GO_BP = function(genes, bg) enrichGO(
    gene = genes, universe = bg,
    OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
    ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE
  ),
  GO_MF = function(genes, bg) enrichGO(
    gene = genes, universe = bg,
    OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
    ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE
  ),
  GO_CC = function(genes, bg) enrichGO(
    gene = genes, universe = bg,
    OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
    ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE
  )
  # Add more  if needed
)
bg_genes <- rownames(dataVoom$E)

all_results <- lapply(names(enrichment_methods), function(method_name) {
  
  method_fun <- enrichment_methods[[method_name]]
  
  res <- lapply(seq_along(genes_by_cluster), function(i) {
    cluster_genes <- genes_by_cluster[[i]]
    method_fun(cluster_genes, bg_genes)
  })
  
  names(res) <- paste0("Cluster_", seq_along(genes_by_cluster))
  res
})


##########
#save all combined

cluster_pr_combined <- lapply(names(all_results), function(db_name) {
  
  cluster_list <- all_results[[db_name]]
  
  lapply(names(cluster_list), function(clust_name) {
    
    x <- cluster_list[[clust_name]]
    
    if (is.null(x)) return(NULL)
    
    df <- tryCatch(x@result, error = function(e) NULL)
    
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df$cluster <- clust_name
    df$db <- db_name
    
    df
  }) %>%
    bind_rows()
  
}) %>%
  bind_rows()
#create directory

write_rds(cluster_pr_combined,paste0(Enrichment_dir,"enrichment_cluster_profiler.rds"))

#per method_percluster--------------
for (method_name in names(all_results)) {
  method_dir <- file.path(Enrichment_dir, method_name)
  if (!dir.exists(method_dir)) dir.create(method_dir)
  
  cluster_list <- all_results[[method_name]]
  
  for (nm in names(cluster_list)) {
    enr <- cluster_list[[nm]]
    
    if (is.null(enr) || !inherits(enr, "enrichResult") || nrow(enr@result) == 0) {
      message("Skipping ", nm, " in ", method_name)
      next
    }
    
    file_out <- file.path(method_dir, paste0(nm, "_dotplot.pdf"))
    
    pdf(file_out, height = 10)
    print(dotplot(enr, showCategory = 15) + ggtitle(paste(nm, method_name)))
    dev.off()
  }
}
#per dbx method
for(dbx in unique(cluster_pr_combined$db)){
  terms <- cluster_pr_combined%>%
    filter(db == dbx)%>%
    filter(p.adjust < 0.05, FoldEnrichment > 15)%>%
    pull(Description)%>%
    unique()
  nrow <- length(terms) + 5
  p <- ggplot(cluster_pr_combined%>%
           filter(Description %in% terms),
         aes(x = gsub("_", " ", cluster), y = Description,
             size = -log10(p.adjust),
             colour = FoldEnrichment))+
    geom_point()+
    scale_color_gradient2(high = "red",
                          mid = "white",
                          low = "blue")+
    labs( x = "gene cluster",
          y = "pathways")+
    optimized_theme_fig()+
    
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12,angle = 45))
    
  ggsave(paste0(Enrichment_dir,"/",dbx,"/enrichment_per_cluster.pdf"), plot = p,
         w = 28, h = min(40,nrow * 2.5), units = "cm" )
}

#enrichR----------------------------------
enrichr_dbs <- c(
  "GO_Biological_Process_2021",
  "KEGG_2021_Mouse",
  "Reactome_2022",
  "MSigDB_Hallmark_2020",
  "WikiPathways_2021_Mouse",
  "TRRUST_Transcription_Factors_2019",
  "GO_Molecular_Function_2023",
  "GO_Biological_Process_2023",
  "CellMarker_2024")


enrichr_results <- lapply(seq_along(genes_by_cluster), function(i) {
  
  cluster_genes <- genes_by_cluster[[i]]
  symbols <- convert_to_symbol(cluster_genes)
  
  if (length(symbols) == 0) return(NULL)
  
  enrichr(symbols, enrichr_dbs)
})
enrichr_long <- map_dfr(seq_along(enrichr_results), function(i) {
  
  cluster_name <- paste0("Cluster ", i)
  res_list <- enrichr_results[[i]]
  
  if (is.null(res_list)) return(NULL)
  
  map_dfr(names(res_list), function(db) {
    
    df <- res_list[[db]]
    
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df %>%
      mutate(
        Cluster = cluster_name,
        Database = db
      )
  })
})

Enrichment_dir <- file.path(cluster_dir, "enrichment_enrichR")
if (!dir.exists(Enrichment_dir)) dir.create(Enrichment_dir)

plot_enrichr <- function(df, title) {
  
  df <- df[order(df$Adjusted.P.value), ]
  df <- head(df, 15)
  
  ggplot(df, aes(x = reorder(Term, -log10(Adjusted.P.value)),
                 y = Cluster,
                 size = -log10(Adjusted.P.value),
                 color = log2(Odds.Ratio))) +
    geom_point() +
    scale_color_gradient(low = "white", high = "red")+
    coord_flip() +
    labs(
      title = title,
      x = "Pathway",
      y = "Clusters",
      size = "-log10 Adjusted P-value",
      color = "log2(Odds.Ratio)"
    )+
    optimized_theme_fig()
}
for (db_name in enrichr_dbs) {
  
  dbs_res <- enrichr_long |> 
    filter(Database == db_name)
  
  
  
  db_dir <- file.path(Enrichment_dir, db_name)
  if (!dir.exists(db_dir)) dir.create(db_dir)
  
  
  
  file_out <- file.path(db_dir, paste0("enrichR", "_", db_name, "_dotplot.pdf"))
  
  p <- plot_enrichr(dbs_res, paste(db_name))
  
  pdf(file_out, height = 8, width = 10)
  print(p)
  dev.off()
}

