#loading libraries-------------------
library(limma)
library(edgeR)
library(tidyverse)
library(ggplot2)
library(readxl)
library(biomaRt)
library(stringr)
library(ggrepel)
library(ComplexHeatmap)
library(patchwork)
library(readxl)
library(dplyr)
library(forcats)
library(circlize)
library(latex2exp)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)  # mouse; use org.Hs.eg.db for human
library(enrichplot)
source("Ag_optimized_theme.R")
# Load data-----------------

data <- read.table("merged_counts.txt", header = TRUE)
metadata <- read_excel("Conditions.xlsx") %>%
  as.data.frame()

# Clean column names in counts data

colnames(data) <- gsub(
  "X.lisc.data.scratch.decker.proseq.proseq_out.|_dedup_QC_end.sort.bam",
  "",
  colnames(data)
)

# Set metadata rownames and ensure Sample column is character

metadata$Sample <- as.character(metadata$Sample)

# Standardize Treatment

metadata <- metadata %>%
  mutate(
    Treatment = gsub("-", "ut", Treatment),
    Treatment = gsub("b", "IFNb", Treatment),
    Treatment = gsub("g", "IFNg", Treatment),
    Treatment = factor(Treatment, levels = c("ut", "IFNb", "IFNg"))
  )

# Standardize Time

metadata <- metadata %>%
  mutate(
    Time = gsub("0", "ut", Time),
    Time = fct_relevel(Time, "ut")
  )

# Relevel Genotype

metadata <- metadata |>
  mutate(
    Genotype = gsub("IRF2-/-","IRF2KO",Genotype))|>
  mutate(Genotype = fct_relevel(Genotype, "WT")
  )


# Create Condition factor

metadata <- metadata %>%
  mutate(
    Condition = paste0(Treatment, "_", Time),
    Condition = fct_relevel(Condition, "ut_ut")
  )

rownames(metadata) <- metadata$Sample
stopifnot(all(rownames(metadata)==colnames(data)))
# correlation plot
corMT <- cor(data)
diag(corMT) <- NA


rownames(corMT)
metadata$samples <- paste0(metadata$Condition, metadata$Genotype, metadata$Replicate)


rownames(corMT) <- metadata$samples[match(rownames(corMT), metadata$Sample)]
colnames(corMT) <- metadata$samples[match(colnames(corMT), metadata$Sample)]
################

#if (!dir.exists("QC_and_basic_plots")) dir.create("QC_and_basic_plots")

#outfile <- file.path("QC_and_basic_plots", "Clustering_of_samples.pdf")
#pdf(outfile, w = 25, h = 25)
Heatmap(corMT, 
        cluster_rows = T, 
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 15))
#dev.off()
#design--------
design <- model.matrix(~Condition*Genotype, data = metadata)

dge <- DGEList(data)
dge <- calcNormFactors(dge, method = "TMM")
keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr,]

#voom
dataVoom <- voom(dge, design=design, plot = TRUE) # insert your model matrix
# PCA
# PCA on voom-transformed data
voom_mat <- dataVoom$E   # log2-CPM values

pca <- prcomp(t(voom_mat), scale. = TRUE)

# Percent variance explained
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1 <- round(percentVar[1], 1)
pc2 <- round(percentVar[2], 1)

# Build PCA dataframe
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Sample = rownames(pca$x)
) %>%
  left_join(metadata, by = c("Sample"))
#
metadata$Condition <- factor(metadata$Condition,
                             levels = c("ut_ut","IFNb_1.5","IFNb_4","IFNb_24", "IFNb_48",
                                        "IFNg_1.5", "IFNg_4", "IFNg_24", "IFNg_48"))
my_colors <- list(
  "grey","#492050", "#82498C", "#B574C2", "#D2A9DB",   
  "#256C26","#91C392", "#4E9D4F","#C8E1C9"
)


# PCA plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, shape = Genotype)) +
  geom_point(size = 4) +
  scale_color_manual(values = my_colors)+
  theme_bw(base_size = 14) +
  labs(
    x = paste0("PC1 (", pc1, "%)"),
    y = paste0("PC2 (", pc2, "%)"),
    title = "PCA on normalized expression data"
  ) 
if (!dir.exists("QC_and_basic_plots")) dir.create("QC_and_basic_plots")

pca <- file.path("QC_and_basic_plots", "PCA_on normalized_expression.pdf")
ggsave(pca)
# -----------------------------
limmaFit <- lmFit(dataVoom, design)
limmaFit <- eBayes(limmaFit)

# Extract results
limmaRes <- map_dfr(colnames(coef(limmaFit)), function(coefx) {
  topTable(limmaFit, coef = coefx, number = Inf) %>%
    rownames_to_column("genes") %>%
    filter(coefx != "(Intercept)") %>%
    mutate(coef = coefx,
           group = case_when(
             logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
             logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
             TRUE ~ "n.s"
           ))
})

# List all treatment conditions excluding baseline (ut_ut)
treatments <- c("IFNb_1.5","IFNb_4","IFNb_24","IFNb_48",
                "IFNg_1.5","IFNg_4","IFNg_24","IFNg_48")

# Step 1: replace ":" with "." in design and limmaFit
colnames(limmaFit$coefficients) <- gsub(":", ".", colnames(limmaFit$coefficients))
colnames(limmaFit$design) <- gsub(":", ".", colnames(limmaFit$design))
if(!is.null(limmaFit$contrasts)) {
  colnames(limmaFit$contrasts) <- gsub(":", ".", colnames(limmaFit$contrasts))
}

# Step 2: generate contrast list using "." instead of ":"
contrast_list <- sapply(treatments, function(trt) {
  main <- paste0("Condition", trt)
  inter <- paste0("Condition", trt, ".GenotypeIRF2KO")
  paste0(main, "+", inter)
})

# Step 3: make contrast matrix using the updated coefficient names
contrast_matrix <- makeContrasts(contrasts = contrast_list,
                                 levels = colnames(limmaFit$coefficients))

# Step 4: fit contrasts
fit2 <- contrasts.fit(limmaFit, contrast_matrix)
fit2 <- eBayes(fit2)

limmaRes_contrasts <- map_dfr(colnames(coef(fit2)), function(coefx) {
  topTable(fit2, coef = coefx, number = Inf) %>%
    rownames_to_column("genes") %>%
    mutate(
      coef = coefx,
      group = case_when(
        logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
        logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
        TRUE ~ "n.s"
      )
    )
})
limmaRes <- bind_rows(limmaRes, limmaRes_contrasts)

limmaRes <- limmaRes %>%
  mutate(
    # Remove prefixes
    coef = str_replace_all(coef, "Condition|Genotype", ""),
    coef = str_replace(coef, "IRF2-/-", "IRF2KO"),
    coef = if_else(
      str_detect(coef, "\\+"),
      str_replace(coef, "\\+.*", "_IRF2KO"), coef ),
    # Handle interactions (contains ":")
    coef = if_else(
      str_detect(coef, ":"),
      # Replace b/g at start with IFNb/IFNg, replace : with _, append _Interaction
      coef %>% str_replace(":", "_") %>%
        paste0("_Interaction"),
      coef
    ))
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ensembl_ids <- unique(limmaRes$genes)
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),  # Ensembl ID → gene symbol
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
limmaRes$gene.name <- gene_info$mgi_symbol[match(limmaRes$genes, gene_info$ensembl_gene_id)]
###
#prepare tables per coef-----------
# Split limmaRes into a list of tibbles, one per coefficient
limmaRes_list <- limmaRes %>%
  group_by(coef) %>%
  group_split() %>%
  setNames(unique(limmaRes$coef))

# Create folder if it doesn't exist
# Split by coefficient (CORRECT & SAFE)
limmaRes_list <- split(limmaRes, limmaRes$coef)

# Create output directory if needed
if (!dir.exists("Table_of_results")) {
  dir.create("Table_of_results")
}

# Save limma results per coefficient
for (coef_name in names(limmaRes_list)) {
  write.csv(
    limmaRes_list[[coef_name]],
    file = file.path(
      "Table_of_results",
      paste0("limma_results_", coef_name, ".csv")
    ),
    row.names = FALSE
  )
}



#Significant results----------------------
limmaRessig <- limmaRes %>%
   filter(adj.P.Val < 0.05, abs(logFC)>1)
# Split SIGNIFICANT results by coefficient
limmaResSig_list <- split(limmaRessig, limmaRessig$coef)

# Create output directory if needed
if (!dir.exists("Table_of_SIGNIFICANT_results")) {
  dir.create("Table_of_SIGNIFICANT_results")
}

# Save significant limma results per coefficient
for (coef_name in names(limmaResSig_list)) {
  write.csv(
    limmaResSig_list[[coef_name]],
    file = file.path(
      "Table_of_SIGNIFICANT_results",
      paste0("limma_SIGNIFICANT_results_", coef_name, ".csv")
    ),
    row.names = FALSE
  )
}

#
p_pval_color <- ggplot(limmaRes, aes(x = P.Value, fill = factor(floor(AveExpr)))) +
  geom_histogram(bins = 50, color = "black", alpha = 0.8) +
  labs(
    title = "P-value distribution colored by average expression",
    x = "Raw p-value",
    fill = "floor(AveExpr)"
  ) +
  facet_wrap(~coef)+
  theme_bw()
if (!dir.exists("QC_and_basic_plots")) dir.create("QC_and_basic_plots")

outfile <- file.path("QC_and_basic_plots", "Pvalue_distribution.pdf")
ggsave(outfile, width = 18, height = 18)
#volcano--------------------

top_genes <- limmaRes %>%
  group_by(coef) %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(desc(abs(logFC))) %>%  # or desc(logFC) if you only care about up
  slice_head(n = 2) %>%
  ungroup()

ggplot() +
  # Hex background for non-significant genes
  stat_bin_hex(
    data = filter(limmaRes, group == "n.s"),
    aes(x = logFC, y = -log10(adj.P.Val), fill = ..count..),
    bins = 20, color = NA, alpha = 0.7
  ) +
  scale_fill_gradient(
    low = "lightgrey",
    high = "black",
    limits = c(1, 5000),
    name = "Gene Count"
  ) +
  
  # Vertical cutoff lines
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.4
  ) +
  
  # Significant up/down genes
  geom_point(
    data = filter(limmaRes, group %in% c("up","down")),
    aes(x = logFC, y = -log10(adj.P.Val), color = group),
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c(
      "up" = "#D0154E",
      "down" = "#4C889C"
    )
  ) +
  
  # Highlight selected genes
  geom_point(
    data = top_genes,
    aes(x = logFC, y = -log10(adj.P.Val)),
    color = "black",
    size = 0.5
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(x = logFC, y = -log10(adj.P.Val), label = gene.name),
    size = 2,
    color = "black",
    max.overlaps = 100,
    force = 10,
    force_pull = 0.1,
    max.iter = 3000,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.color = "black",
    segment.size = 0.3,
    min.segment.length = 0.02,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed", angle = 25)
  ) +
  
  labs(
    title = "DEGs per comparison",
    x = "logFC",
    y = "-log10(adj.P)"
  ) +
  
  facet_wrap(
    ~ factor(coef, levels = c(
      "IFNb_1.5", "IFNb_4", "IFNb_24", "IFNb_48",
      "IFNb_1.5_IRF2KO", "IFNb_4_IRF2KO", "IFNb_24_IRF2KO", "IFNb_48_IRF2KO",
      "IFNb_1.5_IRF2KO_Interaction", "IFNb_4_IRF2KO_Interaction",
      "IFNb_24_IRF2KO_Interaction", "IFNb_48_IRF2KO_Interaction",
      "IFNg_1.5", "IFNg_4", "IFNg_24", "IFNg_48",
      "IFNg_1.5_IRF2KO", "IFNg_4_IRF2KO", "IFNg_24_IRF2KO", "IFNg_48_IRF2KO",
      "IFNg_1.5_IRF2KO_Interaction", "IFNg_4_IRF2KO_Interaction",
      "IFNg_24_IRF2KO_Interaction", "IFNg_48_IRF2KO_Interaction",
      "IRF2KO"
    )),
    ncol = 4,
    scales = "free"
  ) +
  
  optimized_theme_fig()
if(!dir.exists("volcano_plot")) {
  dir.create("volcano_plot")
}
volcanop = file.path("volcano_plot", "volcano_plot_per_condition.pdf")
ggsave(volcanop)
#volcano qc with average expression----------------
ggplot() +
  # Hex background for non-significant genes
  stat_bin_hex(
    data = filter(limmaRes, group == "n.s"),
    aes(x = logFC, y = -log10(adj.P.Val), fill = AveExpr),
    bins = 20, color = NA, alpha = 0.7
  ) +
  scale_fill_gradient(
    low = "lightgrey",
    high = "#1589F0",  # blueish gradient for AveExpr
    name = "AveExpr"
  ) +
  
  # Vertical cutoff lines
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.4
  ) +
  
  # Significant up/down genes colored by AveExpr
  geom_point(
    data = filter(limmaRes, group %in% c("up","down")),
    aes(x = logFC, y = -log10(adj.P.Val), color = AveExpr),
    alpha = 0.8
  ) +
  scale_color_gradient(
    low = "#4C889C",   # low expression → blue
    high = "#D0154E",  # high expression → red
    name = "AveExpr"
  ) +
  labs(
    title = "DEGs per comparison (colored by AveExpr)",
    x = "logFC",
    y = "-log10(adj.P)"
  ) +
  
  facet_wrap(
    ~ factor(coef, levels = c(
      "IFNb_1.5", "IFNb_4", "IFNb_24", "IFNb_48",
      "IFNb_1.5_IRF2KO", "IFNb_4_IRF2KO", "IFNb_24_IRF2KO", "IFNb_48_IRF2KO",
      "IFNb_1.5_IRF2KO_Interaction", "IFNb_4_IRF2KO_Interaction",
      "IFNb_24_IRF2KO_Interaction", "IFNb_48_IRF2KO_Interaction",
      "IFNg_1.5", "IFNg_4", "IFNg_24", "IFNg_48",
      "IFNg_1.5_IRF2KO", "IFNg_4_IRF2KO", "IFNg_24_IRF2KO", "IFNg_48_IRF2KO",
      "IFNg_1.5_IRF2KO_Interaction", "IFNg_4_IRF2KO_Interaction",
      "IFNg_24_IRF2KO_Interaction", "IFNg_48_IRF2KO_Interaction",
      "IRF2KO"
    )),
    ncol = 4,
    scales = "free"
  ) +
  
  optimized_theme_fig()

if(!dir.exists("volcano_plot")) {
  dir.create("volcano_plot")
}
volcanop1 = file.path("volcano_plot", "volcano_plot_per_condition_average_expression.pdf")
ggsave(volcanop1)

#############
selected_coefs <- unique(c(grep("Interaction", limmaRessig$coef, value = TRUE), "IRF2KO"))

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
unique(limmaRessig$coef)


# limmaRessig |>
#   filter(coef == coefx)

#############
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
  if(!dir.exists("logFC_per_comparison")) {
    dir.create("logFC_per_comparison")
  }
  file = file.path("logFC_per_comparison",
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
  if(!dir.exists("Expression_plots_per_comparison")) {
    dir.create("Expression_plots_per_comparison")
  }
  file = file.path("Expression_plots_per_comparison",
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

# ht <- Heatmap(
#   expr_mat_scaled,
#   name = "Expression",
#   col = colorRamp2(
#     c(min(expr_mat_scaled), 0, max(expr_mat_scaled)),
#     c("#4C889C","white","#D0154E")
#   ),
#   cluster_rows = TRUE,
#   cluster_columns = FALSE,
#   clustering_method_rows = "ward.D2",
#   row_km = 4,
#   row_gap = unit(3, "mm"),
#   show_row_names = FALSE,
#   show_column_names = TRUE,
#   top_annotation = col_ha,
#   row_title = "Genes",
#   column_title = "Samples"
# )
library(ComplexHeatmap)
library(circlize)
library(grid)

# -----------------------------
# 1. Set seed for reproducibility
# -----------------------------
set.seed(123)

# -----------------------------
# 2. Define number of clusters
# -----------------------------
number_of_clusters <- 6
cluster_dir <- paste0("Heatmap_clustering", number_of_clusters)

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

# Convert Ensembl IDs to gene names per cluster
gene_names_by_cluster <- lapply(genes_by_cluster, function(ensembl_ids){
  limmaRes$gene.name[match(ensembl_ids, limmaRes$genes)]
})

# Remove any NAs
gene_names_by_cluster <- lapply(gene_names_by_cluster, function(x) x[!is.na(x)])

# Save per-cluster CSVs
if (!dir.exists("Genes_per_cluster")) dir.create("Genes_per_cluster")

for (cl in names(gene_names_by_cluster)) {
  write.csv(
    data.frame(Gene = gene_names_by_cluster[[cl]]),
    file = file.path("Genes_per_cluster", paste0(cl, "_genes.csv")),
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
    scale_color_manual(values = c("WT" = "black", "IRF2KO" = "hotpink4")) +
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

# combined_plot
# if (!dir.exists("Heatmap_clustering")) dir.create("Heatmap_clustering")

# trend <- file.path("Heatmap_clustering",
#                   "Trend_lines.pdf")
ggsave(trend,combined_plot)
######################################
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

# Dotplot for each cluster
dir.create(paste0("Enrichment_plot",cluster_dir), showWarnings = FALSE)

for (nm in names(cluster_enrichment)) {
  
  enr <- cluster_enrichment[[nm]]    # extract the enrichResult object
  
  # Skip if empty or invalid
  if (is.null(enr) || !inherits(enr, "enrichResult")) {
    message("Skipping ", nm, ": not an enrichResult")
    next
  }
  
  file_combined <- file.path(paste0("Enrichment_plot",cluster_dir),
                             paste0(nm, "_Enrichment_plot.pdf"))
  
  pdf(file_combined, h = 10)
  print(dotplot(enr, showCategory = 15) + ggtitle(nm))
  dev.off()
}

ISG_core = read.delim(paste0("Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)


#################################
#saving normalized_table
normalized <- dataVoom$E
colnames(normalized) <- metadata$samples[match(colnames(normalized), metadata$Sample)]
head(normalized)
normalized <- as.data.frame(normalized)
normalized$gene.name <- limmaRes$gene.name[
  match(rownames(normalized), limmaRes$genes)
]

colnames(normalized) <- metadata$samples[match(colnames(normalized), metadata$Sample)]
head(normalized)
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
    scale_color_manual(values = c("WT" = "black", "IRF2KO" = "hotpink4")) +
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

