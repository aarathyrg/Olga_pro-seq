ISG_core = read.delim(paste0("Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)
#####################
ISG_core_id <- limmaRes$genes[match(ISG_core,limmaRes$gene.name)] %>% na.omit()
# --- 2. Subset expression matrix ---

expr_mat <- dataVoom$E[ISG_core_id, ]

# --- 3. Remove rows with all NA or zero variance ---

expr_mat <- expr_mat[rowSums(is.na(expr_mat)) != ncol(expr_mat), ]
expr_mat <- expr_mat[apply(expr_mat, 1, var) != 0, ]
#filtering untreated----------------------------------
ut_cols <- metadata_ordered$Sample[metadata_ordered$Condition == "ut_ut"]

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
# Create output directory
dir.create("ISG_core", showWarnings = FALSE)

# Open PDF device
pdf(file = "ISG_core/Ut.pdf", width = 10, height = 18)

# Create heatmap
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
  row_gap = unit(3, "mm"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  top_annotation = col_ha,
  row_title = "Genes",
  column_title = "Samples"
)

# Draw and save
draw(ht)

# Close device
dev.off()
###################
#If for all samples--------------------
# Use same row orders as for Ut clustering


# --- 2. Subset expression matrix ---

expr_mat <- dataVoom$E[ISG_core_id, ]

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

#rownames
rownames(expr_mat_scaled) <- limmaRes$gene.name[match(rownames(expr_mat_scaled), 
                                                      limmaRes$genes)]
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

# Open PDF device
pdf(file = "ISG_core/All_Conditions.pdf", width = 10, height = 18)

# Create heatmap
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
  row_gap = unit(3, "mm"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  top_annotation = col_ha,
  row_title = "Genes",
  column_title = "Samples"
)

# Draw and save
draw(ht)

# Close device
dev.off()
