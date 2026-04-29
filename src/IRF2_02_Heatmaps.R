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
source("src/Ag_optimized_theme.R")
#############
# Main output directory
outdir <- here("Results/IRF2_02_Heatmaps/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

#load data

metadata <- read_rds(here("Results/IRF2_01_DE/Processed_data/metadata.rds"))
limmaRes <- read_rds(here("Results/IRF2_01_DE/Processed_data/limmaRes.rds"))
dataVoom <- read_rds(here("Results/IRF2_01_DE/Processed_data/dataVoom.rds"))
limmaRessig <- limmaRes %>%
  filter(adj.P.Val < 0.05, abs(logFC)>1)

selected_coefs <- unique(c(grep("Interaction", limmaRessig$coef, value = TRUE),
                           "IRF2KO"))

#function get relevant samples----------
get_relevant_samples <- function(coefx, metadata) {
  
  # Always include baseline samples (ut, any genotype)
  
  baseline_samples <- metadata %>%
    filter(Time == "ut") %>%
    pull(Sample)
  
  if (grepl("_IRF2KO_Interaction", coefx)) {
    # Extract treatment/time prefix from coefficient (before _IRF2KO_Interaction)
    prefix <- sub("_IRF2KO_Interaction", "", coefx)
    
    # Select samples matching this condition prefix EXACTLY OR IRF2KO genotype
    relevant <- metadata %>%
      filter(Condition == prefix) %>%
      pull(Sample)
    
    
  } else if (coefx == "IRF2KO") {
    # For IRF2KO main effect, include all samples
    relevant <- baseline_samples
    
  } else {
    # For other coefficients, include all samples (or adjust if needed)
    relevant <- metadata %>%
      pull(Sample)
  }
  
  # Combine with baseline samples and return unique
  
  unique(c(baseline_samples, relevant))
}
#function get relevant coef-----------
get_related_coefs <- function(coefx) {
  if (coefx == "IRF2KO") return("IRF2KO")
  if (grepl("_IRF2KO_Interaction$", coefx)) {
    base <- sub("_IRF2KO_Interaction$", "", coefx)
    return(c(base, paste0(base, "_IRF2KO"), coefx))
  }
  stop("not valid")
}


#Expression plots_per_coef--------------
for (coefx in selected_coefs) {
  
  # 1. Significant genes
  sig_genes <- limmaRessig$genes[limmaRessig$coef == coefx] %>% unique()
  
  if(length(sig_genes) == 0) next
  
  # 2. Relevant samples
  relevant_samples <- get_relevant_samples(coefx, metadata)
  
  # 3. Expression matrix subset
  expr_df <- scale(dataVoom$E[sig_genes, relevant_samples]) %>% as.data.frame()
  expr_df$gene <- rownames(expr_df)
  
  # 4. Pivot to long format
  expr_long <- expr_df %>%
    pivot_longer(cols = -gene, names_to = "Sample", values_to = "Expression") %>%
    left_join(metadata|>
                filter(Sample %in% relevant_samples), by = c("Sample")) %>%
    mutate(
      gene_name = limmaRessig$gene.name[match(gene, limmaRessig$genes)]
    )
  
  # 5. Define sample order (WT ut → WT trt → IRF2KO ut → IRF2KO trt)
  if (coefx == "IRF2KO") {
    wt_ut <- metadata %>% filter(Genotype == "WT", Time == "ut") %>% pull(samples)
    ko_ut <- metadata %>% filter(Genotype == "IRF2KO", Time == "ut") %>% pull(samples)
    sample_order <- c(wt_ut, ko_ut)
  } else {
    prefix <- sub("_IRF2KO_Interaction", "", coefx)
    wt_ut <- metadata$sample[metadata$Genotype == "WT" & metadata$Time == "ut"]
    ko_ut <- metadata$samples[metadata$Genotype == "IRF2KO" & metadata$Time == "ut"]
    wt_trt <- metadata$samples[metadata$Genotype == "WT" & grepl(paste0("^", prefix), metadata$Condition)]
    ko_trt <- metadata$samples[metadata$Genotype == "IRF2KO" & grepl(paste0("^", prefix), metadata$Condition)]
    
    # wt_ut  <- metadata %>% filter(Genotype == "WT", Time == "ut") %>% pull(samples)
    # wt_trt <- metadata %>% filter(Genotype == "WT", grepl(paste0("^", prefix), Condition)) %>% pull(samples)
    # ko_ut  <- metadata %>% filter(Genotype == "IRF2KO", Time == "ut") %>% pull(samples)
    # ko_trt <- metadata %>% filter(Genotype == "IRF2KO", grepl(paste0("^", prefix), Condition)) %>% pull(samples)
    # 
    sample_order <- unique(c(wt_ut, wt_trt, ko_ut, ko_trt))
  }
  
  unique(expr_long$samples)
  # Now assign factor with correct order
  expr_long$samples <- factor(expr_long$samples, levels = sample_order)
  
  
  expr_long$samples <- factor(expr_long$samples, levels = sample_order)
  
  expr_long <- expr_long |> na.omit()
  
  # ---- Clean expr_long ----
  expr_long <- expr_long |>
    dplyr::filter(!is.na(gene_name))
  
  # ---- ComplexHeatmap heatmap ----
  # ---- heatmap matrix ----
  heat_mat <- expr_long %>%
    dplyr::select(gene_name, samples, Expression) %>%
    pivot_wider(names_from=samples, values_from=Expression) %>%
    tibble::column_to_rownames("gene_name") %>%
    as.matrix()
  
  heat_mat <- heat_mat[, sample_order[sample_order %in% colnames(heat_mat)], drop=FALSE]
  
  col_fun <- circlize::colorRamp2(
    c(min(heat_mat, na.rm=TRUE), 0, max(heat_mat, na.rm=TRUE)),
    c("#4C889C","white","#D0154E")
  )
  
  # --- Check heat_mat
  if (nrow(heat_mat) == 0 || ncol(heat_mat) == 0 || all(is.na(heat_mat))) {
    message("Skipping ", coefx, ": no data for heatmap")
    next
  }
  
  
  # ---- heatmap with FIXED cell size ----
  ht <- Heatmap(
    heat_mat,
    name = "Expression",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_row_dend = F,
    row_names_gp = grid::gpar(fontsize = 7),
    width  = grid::unit(0.7 * ncol(heat_mat), "cm"),  # fixed cell width
    height = grid::unit(0.8 * nrow(heat_mat), "cm"),  # fixed cell height → constant cell size
    column_title = paste("Significant genes:", coefx),
    na_col = "grey"
  )
  ht_grob <- grid.grabExpr(print(ht))
  
  # Wrap into a patchwork‐compatible object
  ht_patch <- wrap_elements(full = ht_grob)
  
  # ---- Determine dynamic height based on number of genes ----
  n_genes <- nrow(heat_mat)
  
  # Minimum and maximum height limits
  base_height <- 5        # for very few genes
  per_gene   <- 0.1       # height per gene
  max_height <- 50        # don't let PDF grow too large
  
  # Compute height
  pdf_height <- max(base_height, min(max_height, base_height + n_genes *0.3))
  
  # ---- extract gene order (without modifying print behavior) ----
  heat_mat <- heat_mat[rowSums(!is.na(heat_mat)) > 0, , drop=FALSE]  # remove all-NA rows
  heat_mat <- heat_mat[apply(heat_mat, 1, var, na.rm=TRUE) != 0, , drop=FALSE]  # remove constant rows
  
  
  ht_drawn <- draw(ht)
  
  # ---- Capture grob AFTER draw() ----
  ht_grob <- grid.grabExpr(draw(ht))
  
  # ---- TRUE gene order from drawn heatmap ----
  ordered_genes <- rownames(heat_mat)[row_order(ht_drawn)]
  
  # ---- DOT PLOT using SAME ORDER ----
  needed <- get_related_coefs(coefx)
  
  dot_df <- limmaRes %>%
    filter(genes %in% sig_genes, coef %in% needed) |>
    mutate(coef = factor(coef, levels = needed))|>
    mutate(gene.name = factor(gene.name, levels = rev(ordered_genes))
    ) |>
    drop_na()
  
  p <- ggplot(dot_df,
              aes(x = coef, y = gene.name,
                  size = pmin(3, -log10(adj.P.Val)),
                  fill = pmax(pmin(logFC, 2), -2))) +  # fill maps to logFC
    geom_point(shape = 21, color = "black") +
    scale_size_continuous(
      range = c(0,3),
      limits = c(0,3),
      #breaks = c(1,3,5),
      name =TeX("$-\\log_{10}(p_{adj})$")
    )+
    scale_fill_gradient2(high="#D0154E", mid="white", low="#4C889C") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(title = paste("Dot plot for", coefx))+
    optimized_theme_fig()+
    theme(axis.text.x = element_text(angle=45, hjust=1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 10, margin = margin(b = 5)),
          plot.margin = margin(t = 22, r = 15, b = 22, l = 15),  # more space around plot
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  
  logFC_plot <- file.path(outdir,"logFC_per_comparison_plots")
  if(!dir.exists(logFC_plot)) {
    dir.create(logFC_plot)
  }
  file = file.path(logFC_plot,
                   paste0(coefx, "dotplot.png"))
  ggsave(filename = file, p, 
         width = pmax(3,length(relevant_samples)*0.5), height=pdf_height, limitsize = F)
  
  
  p1 <- wrap_elements(full = ht_grob) + p
  ht_patch <- wrap_elements(full = ht_grob)
  
  # single row → same height
  
  
  # ---- ensure ggplot respects the same height ----
  #p <- p + coord_fixed(ratio = nrow(heat_mat) / length(sample_order))
  
  # ---- combine side by side with same height ----
  p1 <- ht_patch | p +
    plot_layout(widths = c(1,2), heights = 1)  # side by side, same height
  # Create folder if it doesn't exist
  Expr_plot <-file.path(outdir,"Expression_plots_per_comparison_plots")
  if(!dir.exists(Expr_plot)) {
    dir.create(Expr_plot)
  }
  file = file.path(Expr_plot,
                   paste0(coefx, "heatmap_dotplot.png"))
  ggsave(filename = file, p1, 
         width = length(relevant_samples)*1, height=pdf_height, limitsize = F)
  
}
######################################
#IRF2KO
#Pre-clean Condition into Stim + Time ----
metadata <- metadata %>%
  mutate(
    Stim       = sub("_.*", "", Condition),          # IFNb / IFNg / ut
    Timepoint  = sub(".*_", "", Condition)           # 1.5 / 4 / 24 / 48 / ut
  )

##############################
#heatmap-----------------------
###
#genes of interest


# --- 1. Identify genes of interest (GOI) ---

selected_coefs <- c(grep("Interaction", limmaRessig$coef, value = TRUE), "IRF2KO")
goi <- limmaRessig %>%
  filter(coef %in% selected_coefs) %>%
  pull(genes) %>%
  unique()

# --- 2. Subset expression matrix ---

expr_mat <- dataVoom$E[goi, ]

# --- 3. Remove rows with all NA or zero variance ---

expr_mat <- expr_mat[rowSums(is.na(expr_mat)) != ncol(expr_mat), ]
expr_mat <- expr_mat[apply(expr_mat, 1, var) != 0, ]

# --- 4. Scale expression per gene ---

expr_mat_scaled <- t(scale(t(expr_mat)))

# --- 5. Define color function ---

col_fun <- colorRamp2(
  c(min(expr_mat_scaled), 0, max(expr_mat_scaled)),
  c("#4C889C", "white", "#D0154E")
)


# Make sure columns of expr_mat_scaled match metadata Sample

expr_mat_scaled <- expr_mat_scaled[, colnames(expr_mat_scaled) %in% metadata$Sample]

colnames(expr_mat_scaled) <- metadata$samples[match(metadata$Sample, colnames(expr_mat_scaled))]

# Order metadata to match the heatmap columns
treatments <- c("ut_ut", "IFNb_1.5", "IFNg_1.5", "IFNb_4", "IFNg_4", 
                "IFNb_24", "IFNg_24", "IFNb_48", "IFNg_48")
genotypes <- c("WT1", "WT2", "WT3", "IRF2KO1", "IRF2KO2", "IRF2KO3")

# Use expand.grid to get all combinations in the right order
df <- expand.grid(genotype = genotypes, treatment = treatments, 
                  stringsAsFactors = FALSE)

# Reorder rows to match treatment first, then genotype
df <- df[order(match(df$treatment, treatments), match(df$genotype, genotypes)), ]

# Create sample names
sample_order <- paste0(df$treatment, df$genotype)
# Define colors for Condition


# Reorder the columns of your expression matrix
expr_mat_scaled <- expr_mat_scaled[, sample_order]
# First, create a consistent sample identifier in metadata

# Colorblind-friendly blue-purple palette in the same order as original

condition_colors <- setNames(
  c(
    "darkgrey",
    "#3D1778","#0E3F5C" ,
    "#82498C","#05547F",
    "#B574C2", "#3690C0",
    "#D2A9DB","#6B91C6"
  ),
  unique(metadata$Condition) 
)


ko_colors <- c(
  "WT" = "black",
  "IRF2KO" = "hotpink4"
)
# Optional: reorder your metadata to match
metadata_ordered <- metadata[match(colnames(expr_mat_scaled),metadata$samples),]

# Then create the heatmap
# Column annotation with two layers
metadata_ordered$Genotype <- factor(
  metadata_ordered$Genotype,
  levels = names(ko_colors)   # ensure levels exactly match color names
)

col_ha <- HeatmapAnnotation(
  Genotype = anno_simple(
    as.character(metadata_ordered$Genotype),  # <- convert factor to character
    col = ko_colors,
    height = unit(3, "mm")# thin bar
  ),
  show_annotation_name  = T,
  Condition = anno_simple(
    as.character(metadata_ordered$Condition),
    col = condition_colors,
    height = unit(4, "mm")   # taller bar
  ),
  which = "column",
  show_legend = TRUE
)

# Draw heatmap with annotation
colnames(expr_mat_scaled) <- metadata_ordered$samples
lg_genotype <- Legend(
  title = "Genotype",
  labels = names(ko_colors),
  legend_gp = gpar(fill = ko_colors)
)
lg_condition <- Legend(
  title = "Condition",
  labels = names(condition_colors),
  legend_gp = gpar(fill = condition_colors)
)



# -----------------------------
# 1. Set seed for reproducibility
# -----------------------------
set.seed(123)

# -----------------------------
# 2. Define number of clusters
# -----------------------------
number_of_clusters <- 6
cluster_dir <- paste0(outdir,"Heatmap_clustering", number_of_clusters)

# -----------------------------
# 3. Compute row clusters (k-means)
# -----------------------------
row_clusters <- kmeans(expr_mat_scaled, centers = number_of_clusters, nstart = 50)$cluster
row_clusters <- factor(row_clusters, levels = 1:number_of_clusters)  # ensure consistent ordering

# -----------------------------
# 4. Create Heatmap
# -----------------------------
ht <- Heatmap(
  expr_mat_scaled,
  name = "Expression",
  col = colorRamp2(c(min(expr_mat_scaled), 0, max(expr_mat_scaled)),
                   c("#4C889C", "white", "#D0154E")),
  cluster_rows = TRUE,             # hierarchical clustering within clusters
  clustering_method_rows = "ward.D2",
  row_split = row_clusters,        # fixed clusters
  cluster_columns = FALSE,         # NO column clustering
  row_gap = unit(3, "mm"),
  show_row_names = FALSE,
  show_column_names = TRUE,
  top_annotation = col_ha,
  row_title = "Genes",
  column_title = "Samples"
)

# -----------------------------
# 5. Create output folder
# -----------------------------
if (!dir.exists(cluster_dir)) dir.create(cluster_dir)
heatmap_file <- file.path(cluster_dir, "Interaction_IRF2KO_genes_heatmap_with_condition.pdf")

# -----------------------------
# 6. Save heatmap to PDF
# -----------------------------
pdf(heatmap_file, width = 12, height = 10)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_list = list(lg_genotype, lg_condition),
  annotation_legend_side = "right"
)
dev.off()

cat("Heatmap saved to:", heatmap_file, "\n")


#####cluster gene list-----------

# Draw heatmap and save the object
ht_obj <- draw(ht)
write_rds(ht_obj, paste0(outdir,"/ht_object.rds"))
# Get row indices for each cluster (top to bottom in the heatmap)
row_orders <- row_order(ht_obj)

# Reorder cluster factor levels according to top-to-bottom appearance
cluster_top_to_bottom <- names(row_orders)  # e.g., "1","2",...
row_clusters <- factor(row_clusters, levels = cluster_top_to_bottom)

# Extract gene names for each cluster in top-to-bottom order
genes_by_cluster <- lapply(row_orders, function(idx) {
  rownames(expr_mat_scaled)[idx]
})
names(genes_by_cluster) <- paste0("Cluster", seq_along(genes_by_cluster))
write_rds(genes_by_cluster,paste0(outdir,"/genes_by_cluster.rds"))
# Convert Ensembl IDs to gene names per cluster
gene_names_by_cluster <- lapply(genes_by_cluster, function(ensembl_ids){
  limmaRes$gene.name[match(ensembl_ids, limmaRes$genes)]
})

# Remove any NAs
gene_names_by_cluster <- lapply(gene_names_by_cluster, function(x) x[!is.na(x)])

# Save per-cluster CSVs
Genes_per_cluster <- file.path(cluster_dir,"Genes_per_cluster")
if (!dir.exists(Genes_per_cluster)) dir.create(Genes_per_cluster)

for (cl in names(gene_names_by_cluster)) {
  write.csv(
    data.frame(Gene = gene_names_by_cluster[[cl]]),
    file = file.path(Genes_per_cluster, paste0(cl, "_genes.csv")),
    row.names = FALSE
  )
}

# Check result
gene_names_by_cluster
############
#trend lines-----------------

gene_line_plot <- function(gene_names_by_cluster){
  cluster_name <- names(gene_names_by_cluster)
  genes <- unlist(gene_names_by_cluster)
  
  # Subset and scale
  data <- dataVoom$E[genes,,  drop = FALSE]
  head(data)
  rownames(data) <- limmaRes$gene.name[match(rownames(data), limmaRes$genes)]
  colnames(data) <- metadata$samples[match(colnames(data), metadata$Sample)]
  
  
  # Long format with metadata
  data <- data |>
    as_tibble(rownames = "gene") |>
    pivot_longer(cols = -gene, names_to = "samples", values_to = "E") |>
    left_join(metadata, by = "samples") |>
    group_by(Treatment, Time, Genotype, Condition) |>
    summarise(mean_E = mean(E), .groups = "drop")
  
  # Take only UT rows
  ut <- data |> filter(Treatment == "ut")
  
  # Create two baseline copies
  ut_IFNb <- ut |> mutate(Treatment = "IFNb")
  ut_IFNg <- ut |> mutate(Treatment = "IFNg")
  
  # Combine everything (dropping the original UT rows)
  data2 <- data |> 
    filter(Treatment != "ut") |> 
    bind_rows(ut_IFNb, ut_IFNg) |>
    mutate(
      Treatment = factor(Treatment, levels = c("IFNb", "IFNg")),
      Time = factor(Time, levels = c("ut","1.5","4","24","48"))
    )
  
  
  # ---- Plot ----
  ggplot(data2, aes(x = Time, y = mean_E, color = Genotype, group = Genotype)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("WT" = "black", "IRF2KO" = "darkgreen")) +
    theme_bw() +
    labs(x = "Time (h)",
         y = "Mean expression",
         color = "Genotype",
         title = cluster_name) +
    
    facet_grid(rows = vars(Treatment), scales = "free_y")+optimized_theme_fig()
}
combined_plot <-
  gene_line_plot(genes_by_cluster[1]) +
  gene_line_plot(genes_by_cluster[2]) +
  gene_line_plot(genes_by_cluster[3]) +
  gene_line_plot(genes_by_cluster[4]) +
  gene_line_plot(genes_by_cluster[5]) +
  gene_line_plot(genes_by_cluster[6])

plot_layout(ncol = 2)
trend <- file.path(cluster_dir, "Trend_lines.pdf")


if (!dir.exists(Genes_per_cluster)) dir.create(Genes_per_cluster)

trend <- file.path(Genes_per_cluster,"Trend_lines.pdf")
ggsave(trend,combined_plot)
######################################

ISG_core = read.delim(paste0("Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)



#################################
#saving normalized_table
normalized <- dataVoom$E
colnames(normalized) <- metadata$samples[match(colnames(normalized), metadata$Sample)]

normalized <- as.data.frame(normalized)
normalized$gene.name <- limmaRes$gene.name[
  match(rownames(normalized), limmaRes$genes)
]

colnames(normalized) <- metadata$samples[match(colnames(normalized), metadata$Sample)]

normalized <- as.data.frame(normalized)
normalized$gene.name <- limmaRes$gene.name[
  match(rownames(normalized), limmaRes$genes)
]
head(normalized)
write.csv(normalized, "Normalized_counts.csv")
#example gene plot
# Create a single-gene list
tlr4_list <- list(Tlr4 = "Tlr4")  # the name and gene symbol

# Subset dataVoom$E by Tlr4
gene_line_plot_single <- function(gene_symbol) {
  
  # Map gene symbol to Ensembl ID
  ensembl_id <- limmaRes$genes[match(gene_symbol, limmaRes$gene.name)]
  
  data <- dataVoom$E[ensembl_id,, drop = FALSE]
  rownames(data) <- gene_symbol  # keep symbol for plotting
  colnames(data) <- metadata$Sample[match(colnames(data), metadata$Sample)]
  
  # Long format with metadata
  data_long <- data |>
    as_tibble(rownames = "gene") |>
    pivot_longer(cols = -gene, names_to = "samples", values_to = "E") |>
    left_join(metadata, by = c("samples" = "Sample")) |>
    group_by(Treatment, Time, Genotype, Condition) |>
    summarise(mean_E = mean(E), .groups = "drop")
  
  # Take only UT rows
  ut <- data_long |> filter(Treatment == "ut")
  
  # Create two baseline copies
  ut_IFNb <- ut |> mutate(Treatment = "IFNb")
  ut_IFNg <- ut |> mutate(Treatment = "IFNg")
  
  # Combine everything
  data2 <- data_long |> 
    filter(Treatment != "ut") |> 
    bind_rows(ut_IFNb, ut_IFNg) |>
    mutate(
      Treatment = factor(Treatment, levels = c("IFNb", "IFNg")),
      Time = factor(Time, levels = c("ut","1.5","4","24","48"))
    )
  
  # Plot
  ggplot(data2, aes(x = Time, y = mean_E, color = Genotype, group = Genotype)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("WT" = "black", "IRF2KO" = "darkgreen")) +
    theme_bw() +
    labs(x = "Time (h)",
         y = "Mean expression",
         color = "Genotype",
         title = gene_symbol) +
    facet_grid(rows = vars(Treatment), scales = "free_y") +
    optimized_theme_fig()
}

# Plot Tlr4
tlr4_plot <- gene_line_plot_single("Tlr4")
tlr4_plot
##############
# Count significant genes per coefficient
sig_gene_counts <- limmaRessig %>%
  group_by(coef) %>%
  summarise(
    n_sig_genes = n(),
    .groups = "drop"
  ) 

# View

##################
#enrichment------------------------

# All genes in your dataset
bg_genes <- rownames(dataVoom$E)  # Ensembl IDs

cluster_enrichment <- lapply(1:length(genes_by_cluster), function(i){
  
  cluster_genes <- genes_by_cluster[[i]]  # Ensembl IDs
  
  enrichGO(gene         = cluster_genes,
           universe     = bg_genes,
           OrgDb        = org.Mm.eg.db,
           keyType      = "ENSEMBL",  # Important since you have Ensembl IDs
           ont          = "BP",        # Biological Process
           pAdjustMethod= "BH",
           qvalueCutoff = 0.05,
           readable     = TRUE)        # converts Ensembl IDs to gene symbols
})

names(cluster_enrichment) <- paste0("Cluster_", 1:length(genes_by_cluster))



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
  ),
  KEGG = function(genes, bg) enrichKEGG(
    gene = genes,
    universe = bg,
    organism = "mmu",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  # Add more (e.g. Reactome) if needed
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

names(all_results) <- names(enrichment_methods)
Enrichment_dir <- file.path(cluster_dir, "enrichment_cluster_profiler")
if (!dir.exists(Enrichment_dir)) dir.create(Enrichment_dir)

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
# Dotplot for each cluster
Enrichment_dir <- file.path(cluster_dir,"enrichment_cluster_profiler/")
if (!dir.exists(Enrichment_dir)) dir.create(Enrichment_dir)


for (nm in names(cluster_enrichment)) {
  
  enr <- cluster_enrichment[[nm]]    # extract the enrichResult object
  
  # Skip if empty or invalid
  if (is.null(enr) || !inherits(enr, "enrichResult")) {
    message("Skipping ", nm, ": not an enrichResult")
    next
  }
  
  file_combined <- file.path(Enrichment_dir,
                             paste0(nm, "_Enrichment_plot.pdf"))
  
  pdf(file_combined, h = 10)
  print(dotplot(enr, showCategory = 15) + ggtitle(nm))
  dev.off()
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
  
  if (is.null(cluster_res)) next
  
  db_dir <- file.path(Enrichment_dir, db_name)
  if (!dir.exists(db_dir)) dir.create(db_dir)
  
  
    
  file_out <- file.path(db_dir, paste0("enrichR", "_", db_name, "_dotplot.pdf"))
    
  p <- plot_enrichr(dbs_res, paste(db_name))
    
    pdf(file_out, height = 8, width = 10)
    print(p)
    dev.off()
}

