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

#limmaRes
limmaRes <- read_rds("limmaRes_IRF2.rds")
unique(limmaRes$coef)
head(limmaRes)
#sinificant
genes_ut <- limmaRes |>
  filter(coef == "IRF2KO")|>
  filter(group != "n.s") |>
  pull(gene.name)
ISG_core = read.delim(paste0("Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)


selected_list <- intersect(genes_ut,ISG_core)
#selected_list <- genes_ut

#####################
selected_list_id <- limmaRes$genes[match(selected_list,limmaRes$gene.name)] %>% na.omit()
# --- 2. Subset expression matrix ---

expr_mat <- dataVoom$E[selected_list_id, ]

# --- 3. Remove rows with all NA or zero variance ---

expr_mat <- expr_mat[rowSums(is.na(expr_mat)) != ncol(expr_mat), ]
expr_mat <- expr_mat[apply(expr_mat, 1, var) != 0, ]
#filtering untreated----------------------------------
ut_cols <- metadata$Sample[metadata$Condition == "ut_ut"]

expr_scaled_ut <- expr_mat[, ut_cols]

# --- 4. Scale expression per gene ---
expr_scaled_ut <- t(scale(t(expr_scaled_ut)))
expr_mat_scaled <- expr_scaled_ut

# --- 5. Define color function ---

col_fun <- colorRamp2(
  c(min(expr_mat_scaled), 0, max(expr_mat_scaled)),
  c("#4C889C", "white", "#D0154E")
)


# Make sure columns of expr_mat_scaled match metadata Sample

expr_mat_scaled <- expr_mat_scaled[, colnames(expr_mat_scaled) %in% metadata$Sample]

colnames(expr_mat_scaled) <- metadata$samples[match(colnames(expr_mat_scaled),metadata$Sample)]
# Order metadata to match the heatmap columns
treatments <- c("ut_ut")
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
rownames(expr_mat_scaled) <- limmaRes$gene.name[match(rownames(expr_mat_scaled), 
                                                      limmaRes$genes)]
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

ht <- Heatmap(
  expr_mat_scaled,
  name = "Expression",
  col = colorRamp2(
    c(min(expr_mat_scaled), 0, max(expr_mat_scaled)),
    c("#4C889C","white","#D0154E")
  ),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_method_rows = "ward.D2",
  #row_km = 4,
  row_gap = unit(3, "mm"),
  show_row_names = T,
  show_column_names = TRUE,
  top_annotation = col_ha,
  row_title = "Genes",
  column_title = "Samples"
)


# Suppose your expression matrix is 'expr_mat'
# Column names should indicate sample/replicate, e.g., "WT_1", "WT_2", "IRF2KO_1", "IRF2KO_2"
colnames(dataVoom$E)
expr <- dataVoom$E
colnames(expr) <- metadata$samples[match(colnames(expr),metadata$Sample)]
sample_IRF2KO2 <- "ut_utIRF2KO2"
sample_WT2 <- "ut_utWT2"

# Compute difference from mean of other samples
diff_IRF2KO2 <- expr[, sample_IRF2KO2] - rowMeans(expr[, colnames(expr) != sample_IRF2KO2])
diff_WT2     <- expr[, sample_WT2] - rowMeans(expr[, colnames(expr) != sample_WT2])

# Quick summary
summary(diff_IRF2KO2)
summary(diff_WT2)

# Optional: histogram to visualize
hist(diff_IRF2KO2, breaks=50, main="IRF2KO Rep 2 vs others", xlab="Difference")
hist(diff_WT2, breaks=50, main="WT Rep 2 vs others", xlab="Difference")

cor_IRF2KO2 <- cor(expr[, sample_IRF2KO2], expr[, grep("IRF2KO", colnames(expr))[-2]])
cor_WT2     <- cor(expr[, sample_WT2], expr[, grep("WT", colnames(expr))[-2]])

cor_IRF2KO2
cor_WT2
