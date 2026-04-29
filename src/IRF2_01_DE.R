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
library(readxl)
library(dplyr)
library(forcats)
library(circlize)
library(latex2exp)
library(patchwork)
library(clusterProfiler)
library(renv)
library(org.Mm.eg.db)  # mouse; use org.Hs.eg.db for human
library(enrichplot)
library(here)
source("src/Ag_optimized_theme.R")
# Load data-----------------

# Main output directory
outdir <- here("Results/IRF2_01_DE")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

#load data
data <- read.table(here("rawdata", "merged_counts.txt"), header = TRUE)

metadata <- read_excel(here("rawdata", "Conditions.xlsx")) %>%
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
pr_data <- file.path(outdir, "Processed_data")
if (!dir.exists(pr_data)) dir.create(pr_data, recursive = TRUE)
outfile <- file.path(pr_data, "metadata.rds")
metadata %>% write_rds(outfile)

rownames(corMT) <- metadata$samples[match(rownames(corMT), metadata$Sample)]
colnames(corMT) <- metadata$samples[match(colnames(corMT), metadata$Sample)]
################
qcdir <- file.path(outdir, "QC_and_basic_plots")
if (!dir.exists(qcdir)) dir.create(qcdir, recursive = TRUE)

outfile <- file.path(qcdir, "Clustering_of_samples.pdf")
pdf(outfile, w = 25, h = 25)
Heatmap(corMT, 
        cluster_rows = T, 
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 15))
dev.off()
#design--------
design <- model.matrix(~Condition*Genotype, data = metadata)

dge <- DGEList(data)
dge <- calcNormFactors(dge, method = "TMM")
keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr,]

#voom
dataVoom <- voom(dge, design=design, plot = TRUE) # insert your model matrix
#save
pr_data <- file.path(outdir, "Processed_data")
if (!dir.exists(pr_data)) dir.create(pr_data, recursive = TRUE)

outfile <- file.path(pr_data, "dataVoom.rds")

dataVoom %>% write_rds(outfile)
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
qcdir <- file.path(outdir, "QC_and_basic_plots")
if (!dir.exists(qcdir)) dir.create(qcdir, recursive = TRUE)
pca <- file.path(qcdir, "PCA_on normalized_expression.pdf")
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
pr_data <- file.path(outdir, "Processed_data")
if (!dir.exists(pr_data)) dir.create(pr_data, recursive = TRUE)

limma_file <- file.path(pr_data, "limmaRes.rds")
write_rds(limmaRes,limma_file)
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
res_dir <- file.path(pr_data, "Table_of_results")


if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}

# Save limma results per coefficient
for (coef_name in names(limmaRes_list)) {
  write.csv(
    limmaRes_list[[coef_name]],
    file = file.path(
      res_dir,
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
sig_res_dir <- file.path(pr_data, "Table_of_SIGNIFICANT_results")


if (!dir.exists(sig_res_dir)) {
  dir.create(sig_res_dir)
}

# Save significant limma results per coefficient
for (coef_name in names(limmaResSig_list)) {
  write.csv(
    limmaResSig_list[[coef_name]],
    file = file.path(
      sig_res_dir,
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
qcdir <- file.path(outdir, "QC_and_basic_plots")
if (!dir.exists(qcdir)) dir.create(qcdir, recursive = TRUE)

outfile <- file.path(qcdir, "P_val_dist.pdf")
ggsave(outfile,plot = p_pval_color, width = 18, height = 18)

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
volcanodir <- file.path(outdir, "Volcano_plots")
if (!dir.exists(volcanodir)) dir.create(volcanodir, recursive = TRUE)

volcanop = file.path(volcanodir, "volcano_plot_per_condition.pdf")
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


volcanop1 = file.path(volcanodir, "volcano_plot_per_condition_average_expression.pdf")
ggsave(volcanop1)

