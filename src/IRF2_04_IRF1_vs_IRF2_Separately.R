################################################################################
## LIBRARIES
################################################################################
library(tidyverse)
library(readxl)
library(forcats)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(here)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
source("Ag_optimized_theme.R")
################################################################################
outdir <- here("Results/IRF2_05_IRF1_vs_IRF2_Separately/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

pr_data <- file.path(outdir, "Processed_data")
if (!dir.exists(pr_data)) dir.create(pr_data, recursive = TRUE)
heatmap_both <- file.path(outdir, "Heatmap_both")
if (!dir.exists(heatmap_both)) dir.create(heatmap_both, recursive = TRUE)
##############################################################################
IRF1_limmaRes <- read_rds(here("Results/IRF2_04_reanalyze_IRF1/Processed_data/limmaRes_IRF1.rds"))
IRF1_metadata <- read_rds(here("Results/IRF2_04_reanalyze_IRF1/Processed_data/metadata.rds"))
IRF1_dataVoom <- read_rds(here("Results/IRF2_04_reanalyze_IRF1/Processed_data/dataVoom.rds"))
IRF2_metadata <- read_rds(here("Results/IRF2_01_DE/Processed_data/metadata.rds"))
IRF2_limmaRes <- read_rds(here("Results/IRF2_01_DE/Processed_data/limmaRes.rds"))
IRF2_dataVoom <- read_rds(here("Results/IRF2_01_DE/Processed_data/dataVoom.rds"))

#
colnames(IRF2_dataVoom$E) <- IRF2_metadata$samples[
  match(colnames(IRF2_dataVoom$E),IRF2_metadata$Sample)
]
rownames(IRF2_dataVoom$E) <- IRF2_limmaRes$gene.name[
  match(rownames(IRF2_dataVoom$E), IRF2_limmaRes$genes)
]
IRF2_metadata
####################################
WT <- IRF2_limmaRes

WT <- WT |>
  filter(!(grepl("IRF2",WT$coef)))


WT <- WT |>
  mutate(Condition = coef)|>
  dplyr::select(
    gene.name,
    Condition,
    logFC_WT = logFC,
    P_WT     = P.Value,
    adjP_WT  = adj.P.Val
  )
unique(WT$Condition)
## 1. Split within-dataset interaction results
unique(IRF2_limmaRes$coef)
unique(IRF1_limmaRes$coef)
res_IRF2 <- IRF2_limmaRes|>
  filter(
    coef == "IRF2KO" |
      grepl("_IRF2KO_Interaction$", coef)
  )|>
  mutate(
    Condition = ifelse(
      coef == "IRF2KO",
      "ut",
      gsub("_IRF2KO_Interaction$", "", coef)
    )
  )



res_IRF1 <- IRF1_limmaRes|>
  filter(
    coef == "IRF1KO" |
      grepl("_IRF1KO_Interaction$", coef)
  )|>
  mutate(
    Condition = ifelse(
      coef == "IRF1KO",
      "ut",
      gsub("_IRF1KO_Interaction$", "", coef)
    )
  )
res_IRF1$gene.name <- res_IRF1$genes
#save_ut for later
res_IRF2_ut <- res_IRF2 |> filter(Condition == "ut")
res_IRF1_ut <- res_IRF1 |> filter(Condition == "ut")
ut <- res_IRF2_ut|>
  dplyr::select(
    gene.name,
    Condition,
    logFC_IRF2KO = logFC,
    P_IRF2     = P.Value,
    adjP_IRF2KO  = adj.P.Val
  )|>
  inner_join(
    res_IRF1_ut|>
      dplyr::select(
        gene.name,
        Condition,
        logFC_IRF1KO = logFC,
        P_IRF1     = P.Value,
        adjP_IRF1KO  = adj.P.Val
      ),
    by = c("gene.name", "Condition")
  )|>
  mutate(group = ifelse(abs(logFC_IRF2KO) > 0.5 & 
                          abs(logFC_IRF1KO) > 0.5 &
                          
                          adjP_IRF2KO < 0.05 &
                          adjP_IRF1KO < 0.05, "both",
                        ifelse(abs(logFC_IRF1KO) > 0.5 & adjP_IRF1KO < 0.05 ,"IRF1", 
                               ifelse(abs(logFC_IRF2KO) > 0.5 & adjP_IRF2KO < 0.05 ,"IRF2" ,"ns"))))

head(ut)
##############
#without WT filtering-
res_joined_MOD <- res_IRF2|>
  dplyr::select(
    gene.name,
    Condition,
    logFC_IRF2KO = logFC,
    P_IRF2     = P.Value,
    adjP_IRF2KO  = adj.P.Val
  )|>
  inner_join(
    res_IRF1|>
      dplyr::select(
        gene.name,
        Condition,
        logFC_IRF1KO = logFC,
        P_IRF1     = P.Value,
        adjP_IRF1KO  = adj.P.Val
      ),
    by = c("gene.name", "Condition")
  )|>
  mutate(group = ifelse(abs(logFC_IRF2KO) > 0.5 & 
                          abs(logFC_IRF1KO) > 0.5 &
                          adjP_IRF2KO < 0.05 &
                          adjP_IRF1KO < 0.05, "both",
                        ifelse(abs(logFC_IRF1KO) > 0.5 & adjP_IRF1KO < 0.05,"IRF1", 
                               ifelse(abs(logFC_IRF2KO) > 0.5 & adjP_IRF2KO < 0.05 
                                      ,"IRF2" ,"ns"))))

#add ut----------
res_joined_MOD <- res_joined_MOD |>
  #dplyr::select(-c("logFC_WT","P_WT","adjP_WT"))|>
  rbind(ut)


res_joined_MOD$Condition <- factor(res_joined_MOD$Condition, levels = c("ut", "IFNb_1.5","IFNb_4", "IFNb_24", "IFNb_48",
                                                                        "IFNg_1.5", "IFNg_4", "IFNg_24", "IFNg_48"))

unique(res_joined_MOD$Condition)
write_rds(res_joined_MOD,paste0(pr_data,"/res_joined_IRF1_IRF2_significant_genes_WITHOUT_WT.rds"))
ggplot(res_joined_MOD, aes(x = logFC_IRF2KO, y = logFC_IRF1KO)) +
  
  # 1️⃣ Background hex for ns
  geom_hex(
    data = dplyr::filter(res_joined_MOD, group == "ns"),
    bins = 50,
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "grey", high = "grey40", guide = "none") +
  
  # 2️⃣ IRF1 & IRF2 first
  geom_point(
    data = dplyr::filter(res_joined_MOD, group %in% c("IRF1", "IRF2")),
    aes(color = group),
    alpha = 0.7,
    size = 1.2
  ) +
  
  # 3️⃣ BOTH last (strictly on top)
  geom_point(
    data = dplyr::filter(res_joined_MOD, group == "both"),
    aes(color = group),
    alpha = 0.9,
    size = 1.8
  ) +
  
  scale_color_manual(
    values = c(
      IRF2 = "#1C3C13",
      IRF1 = "#4A3563",
      both = "#A24C06"
    ),
    name = "Significant in"
  ) +
  
  coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
  facet_wrap(~Condition) +
  optimized_theme_fig() +
  theme(panel.grid.major = NULL)
#save
volcano_dir <- file.path(outdir, "Comparison_volcano_plots")
if (!dir.exists(volcano_dir)) dir.create(volcano_dir, recursive = TRUE)
ggsave(paste0(volcano_dir,"/significant_genes_IRF1_IRF2_both_WITHOUT_WT_induction.pdf"))

## 2. Join IRF2KO (d1) and IRF1KO (d2) effects per gene + condition
res_joined_MOD <- res_IRF2|>
  dplyr::select(
    gene.name,
    Condition,
    logFC_IRF2KO = logFC,
    P_IRF2     = P.Value,
    adjP_IRF2KO  = adj.P.Val
  )|>
  inner_join(
    res_IRF1|>
      dplyr::select(
        gene.name,
        Condition,
        logFC_IRF1KO = logFC,
        P_IRF1     = P.Value,
        adjP_IRF1KO  = adj.P.Val
      ),
    by = c("gene.name", "Condition")
  )%>%
  inner_join(
    WT|>
      dplyr::select(
        gene.name,
        Condition,
        logFC_WT,
        P_WT,
        adjP_WT
      ),
    by = c("gene.name", "Condition"))|>
  mutate(group = ifelse(abs(logFC_IRF2KO) > 0.5 & 
                          abs(logFC_IRF1KO) > 0.5 &
                          abs(logFC_WT) >0.5 &
                          adjP_WT < 0.05 &
                          adjP_IRF2KO < 0.05 &
                          adjP_IRF1KO < 0.05, "both",
                        ifelse(abs(logFC_IRF1KO) > 0.5 & adjP_IRF1KO < 0.05 & adjP_WT < 0.05 & abs(logFC_WT)> 0.5,"IRF1", 
                               ifelse(abs(logFC_IRF2KO) > 0.5 & adjP_IRF2KO < 0.05 & adjP_WT < 0.05 & abs(logFC_WT)> 0.5,"IRF2" ,"ns"))))


unique(res_joined_MOD$Condition)
head(res_joined_MOD)

#add ut----------
res_joined_MOD <- res_joined_MOD |>
  dplyr::select(-c("logFC_WT","P_WT","adjP_WT"))|>
  rbind(ut)
  

res_joined_MOD$Condition <- factor(res_joined_MOD$Condition, levels = c("ut", "IFNb_1.5","IFNb_4", "IFNb_24", "IFNb_48",
                                                                        "IFNg_1.5", "IFNg_4", "IFNg_24", "IFNg_48"))

unique(res_joined_MOD$Condition)
write_rds(res_joined_MOD,paste0(pr_data,"/res_joined_IRF1_IRF2_significant_genes_WITH_WT_induction.rds"))
ggplot(res_joined_MOD, aes(x = logFC_IRF2KO, y = logFC_IRF1KO)) +
  
  # 1️⃣ Background hex for ns
  geom_hex(
    data = dplyr::filter(res_joined_MOD, group == "ns"),
    bins = 50,
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "grey", high = "grey40", guide = "none") +
  
  # 2️⃣ IRF1 & IRF2 first
  geom_point(
    data = dplyr::filter(res_joined_MOD, group %in% c("IRF1", "IRF2")),
    aes(color = group),
    alpha = 0.7,
    size = 1.2
  ) +
  
  # 3️⃣ BOTH last (strictly on top)
  geom_point(
    data = dplyr::filter(res_joined_MOD, group == "both"),
    aes(color = group),
    alpha = 0.9,
    size = 1.8
  ) +
  
  scale_color_manual(
    values = c(
      IRF2 = "#1C3C13",
      IRF1 = "#4A3563",
      both = "#A24C06"
    ),
    name = "Significant in"
  ) +
  
  coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
  facet_wrap(~Condition) +
  optimized_theme_fig() +
  theme(panel.grid.major = NULL)


#save
volcano_dir <- file.path(outdir, "Comparison_volcano_plots")
if (!dir.exists(volcano_dir)) dir.create(volcano_dir, recursive = TRUE)
ggsave(paste0(volcano_dir,"/significant_genes_IRF1_IRF2_both_considering_WT_induction.pdf"))

###################################
#per group enrichment------------------
# Remove ns if you don’t want it
res_both <- res_joined_MOD|> filter(group != "ns")%>%
  filter(group == "both")

#function to plot--------------
plot_condition_heatmap <- function(cond_name,
                                   res_both,
                                   IRF1_dataVoom,
                                   IRF2_dataVoom,
                                   IRF1_metadata,
                                   IRF2_metadata,
                                   outdir) {
  
  message("Processing: ", cond_name)
  
  # ---- 1. Genes significant in this condition ----
  genes_cond <- res_both |>
    dplyr::filter(Condition == cond_name) |>
    dplyr::pull(gene.name) |>
    unique()
  
  if (length(genes_cond) == 0) {
    message("No genes for ", cond_name)
    return(NULL)
  }
  
  # ---- 2. Subset metadata ----
  meta_IRF1 <- IRF1_metadata |>
    dplyr::filter(Condition %in% c("ut_ut", cond_name)) |>
    dplyr::mutate(
      Condition = factor(Condition, levels = c("ut_ut", cond_name)),
      Genotype  = factor(Genotype, levels = c("WT", "IRF1KO"))
    ) |>
    dplyr::arrange(Genotype, Condition, Replicate)
  
  meta_IRF2 <- IRF2_metadata |>
    dplyr::filter(Condition %in% c("ut_ut", cond_name)) |>
    dplyr::mutate(
      Condition = factor(Condition, levels = c("ut_ut", cond_name)),
      Genotype  = factor(Genotype, levels = c("WT", "IRF2KO"))
    ) |>
    dplyr::arrange(Genotype, Condition, Replicate)
  
  samples_IRF1 <- meta_IRF1$samples
  samples_IRF2 <- meta_IRF2$samples
  
  # ---- 3. Expression matrices ----
  # ---- safe gene/sample intersection ----
  
  expr_IRF1 <- IRF1_dataVoom$E
  expr_IRF2 <- IRF2_dataVoom$E
  
  mat_IRF1 <- expr_IRF1[genes_cond, samples_IRF1, drop = FALSE]
  mat_IRF2 <- expr_IRF2[genes_cond, samples_IRF2, drop = FALSE]
  
  # ---- 4. Clean missing genes ----
  mat_IRF1 <- mat_IRF1[complete.cases(mat_IRF1), , drop = FALSE]
  mat_IRF2 <- mat_IRF2[complete.cases(mat_IRF2), , drop = FALSE]
  
  if (nrow(mat_IRF1) == 0 || nrow(mat_IRF2) == 0) {
    message("Empty matrix for ", cond_name)
    return(NULL)
  }
  
  # ---- 5. Combine IRF1 + IRF2 side-by-side ----
  # ---- 5. Z-score within each dataset separately ----
  mat_IRF1_z <- t(scale(t(mat_IRF1)))
  mat_IRF2_z <- t(scale(t(mat_IRF2)))
  
  mat_IRF1_z[is.na(mat_IRF1_z)] <- 0
  mat_IRF2_z[is.na(mat_IRF2_z)] <- 0
  
  # ---- 6. Combine AFTER scaling ----
  mat_z <- cbind(mat_IRF1_z, mat_IRF2_z)
  # ---- 7. Annotation ----
  meta_all <- rbind(
    meta_IRF1 |> dplyr::mutate(Dataset = "IRF1"),
    meta_IRF2 |> dplyr::mutate(Dataset = "IRF2")
  )
  
  cond_colors <- c("ut_ut" = "#ABC5BB")
  cond_colors[cond_name] <- "#3C549A"
  
  col_ha <- ComplexHeatmap::HeatmapAnnotation(
    Condition = meta_all$Condition,
    Genotype  = meta_all$Genotype,
    Dataset   = meta_all$Dataset,
    col = list(
      Condition = cond_colors,
      Genotype = c(WT = "grey60",
                   IRF1KO = "#4A3563",
                   IRF2KO = "darkgreen"),
      Dataset = c(IRF1 = "#E69F00", IRF2 = "#56B4E9")
    )
  )
  
  # ---- 7. FIXED CELL SIZE ----
  cell_w <- unit(4, "mm")
  cell_h <- unit(4, "mm")
  
  ht_width  <- cell_w * ncol(mat_z)
  ht_height <- cell_h * nrow(mat_z)
  
  ht <- ComplexHeatmap::Heatmap(
    mat_z,
    name = "z-score",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    top_annotation = col_ha,
    column_split = meta_all$Dataset,
    show_row_names = TRUE,
    show_column_names = FALSE,
    column_title = paste("BOTH IRFs —", cond_name),
    width = ht_width,
    height = ht_height,
    col = circlize::colorRamp2(
      c(-2, 0, 2),
      c("#4C889C", "white", "#D0154E")
    )
  )
  
  # ---- 8. OUTPUT ----
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  pdf(file.path(outdir, paste0("Heatmap_both_", cond_name, ".pdf")),
      width = 8,
      height = 10)
  
  ComplexHeatmap::draw(ht)
  dev.off()

}
conditions_to_plot <- c(
  "IFNb_1.5", "IFNb_4",
  "IFNg_4", "IFNg_24", "IFNg_48"
)

lapply(conditions_to_plot, function(cn) {
  plot_condition_heatmap(
    cond_name = cn,
    res_both = res_both,
    IRF1_dataVoom = IRF1_dataVoom,
    IRF2_dataVoom = IRF2_dataVoom,
    IRF1_metadata = IRF1_metadata,
    IRF2_metadata = IRF2_metadata,
    outdir = heatmap_both
  )
})
