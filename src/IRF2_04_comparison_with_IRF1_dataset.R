################################################################################
## LIBRARIES
################################################################################
library(tidyverse)
library(limma)
library(edgeR)
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
outdir <- here("Results/IRF2_04_comparison_with_IRF1_dataset/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
################################################################################
## DATASET 1
################################################################################

## ---- Load counts + metadata ----
data1 <- read.table(here("rawdata/IRF1_IRF2/","merged_counts.txt"), header = TRUE)

metadata1 <- read_excel(here("rawdata/IRF1_IRF2/","Conditions.xlsx")) %>%
  as.data.frame()

## ---- Clean count column names ----
colnames(data1) <- gsub(
  "X.lisc.data.scratch.decker.proseq.proseq_out.|_dedup_QC_end.sort.bam",
  "",
  colnames(data1)
)

metadata1$Sample <- as.character(metadata1$Sample)

## ---- Standardize metadata ----
metadata1 <- metadata1 %>%
  mutate(
    Treatment = Treatment %>%
      gsub("-", "ut", .) %>%
      gsub("b", "IFNb", .) %>%
      gsub("g", "IFNg", .),
    Treatment = factor(Treatment, levels = c("ut", "IFNb", "IFNg")),

    Time = gsub("0", "ut", Time),
    Time = fct_relevel(Time, "ut"),

    Genotype = gsub("IRF2-/-", "IRF2KO", Genotype),
    Genotype = fct_relevel(Genotype, "WT"),

    Condition = paste0(Treatment, "_", Time),
    Condition = fct_relevel(Condition, "ut_ut"),

    samples = paste0(Condition, Genotype, Replicate)
  )

## ---- Enforce contracts ----
rownames(metadata1) <- metadata1$Sample
stopifnot(all(colnames(data1) == metadata1$Sample))

colnames(data1) <- metadata1$samples
rownames(metadata1) <- metadata1$samples

################################################################################
## DATASET 2
################################################################################

## ---- Load counts ----
data2 <- read.table(here("rawdata/IRF1_IRF2/",
                         "counts_proseq_correct_index_11_07_22_new.txt"), 
                    header = TRUE)

## ---- Keep WT + IRF1 only ----
keep_cols <- grepl("^(WT|IRF1)_", colnames(data2))
data2 <- data2[, keep_cols]

## ---- Build metadata from column names ----
metadata2 <- tibble(Sample = colnames(data2))

parts <- strsplit(metadata2$Sample, "_")

metadata2 <- metadata2 %>%
  mutate(
    Genotype  = sapply(parts, `[`, 1),
    Time      = sapply(parts, `[`, 2),
    Treatment = sapply(parts, `[`, 3),
    Replicate = sapply(parts, `[`, 4)
  ) %>%
  mutate(
    Genotype = ifelse(Genotype == "IRF1", "IRF1KO", "WT"),
    Genotype = factor(Genotype, levels = c("WT", "IRF1KO")),

    Replicate = gsub("^R", "", Replicate),

    Time = gsub("1h30min", "1.5", Time),
    Time = gsub("h", "", Time),
    Time = factor(Time, levels = c("ut", "1.5", "4", "24", "48")),

    Treatment = factor(Treatment, levels = c("ut", "IFNb", "IFNg")),

    Condition = paste0(Treatment, "_", Time),
    Condition = fct_relevel(Condition, "ut_ut"),

    samples = paste0(Condition, Genotype, Replicate)
  )

## ---- Enforce contracts ----
metadata2 <- as.data.frame(metadata2)
rownames(metadata2) <- metadata2$samples

colnames(data2) <- metadata2$samples[
  match(colnames(data2), metadata2$Sample)
]

rownames(metadata2) <- metadata2$samples

stopifnot(all(colnames(data2) == metadata2$samples))

################################################################################
## COMBINE DATASETS
################################################################################

## ---- Align genes ----
common_genes <- intersect(rownames(data1), rownames(data2))

data_all <- cbind(
  data1[common_genes, ],
  data2[common_genes, ]
)
metadata1 <- metadata1 %>%
  mutate(Replicate = as.character(Replicate))

metadata2 <- metadata2 %>%
  mutate(Replicate = as.character(Replicate))
metadata1 <- metadata1 %>%
  mutate(Dataset = "d1")

metadata2 <- metadata2 %>%
  mutate(Dataset = "d2")
metadata1 <- metadata1 %>%
  mutate(samples = paste0(Dataset, "_", samples))
metadata2 <- metadata2 %>%
  mutate(samples = paste0(Dataset, "_", samples))

colnames(data2) <- metadata2$samples
rownames(metadata2) <- metadata2$samples

colnames(data1) <- metadata1$samples
rownames(metadata1) <- metadata1$samples
stopifnot(
  !any(duplicated(metadata1$samples)),
  !any(duplicated(metadata2$samples))
)

metadata_all <- bind_rows(metadata1, metadata2)
rownames(metadata_all) <- metadata_all$samples
metadata_all <- bind_rows(metadata1, metadata2)
rownames(metadata_all) <- metadata_all$samples
########

ensembl_ids <- rownames(data1)

gene_map <- mapIds(
  org.Mm.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
data1$GeneSymbol <- gene_map
data1 <- data1[!is.na(data1$GeneSymbol), ]
data1 <- data1 %>%
  as.data.frame() %>%
  rowwise() %>%
  mutate(mean_expr = mean(c_across(where(is.numeric)))) %>%
  ungroup() %>%
  arrange(desc(mean_expr)) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  dplyr::select(-mean_expr)
data1 <- as.data.frame(data1)
rownames(data1) <- data1$GeneSymbol
data1$GeneSymbol <- NULL

common_genes <- intersect(rownames(data1), rownames(data2))

data_all <- cbind(
  data1[common_genes, ],
  data2[common_genes, ]
)

stopifnot(all(colnames(data_all) == metadata_all$samples))

## ---- Final checks ----
stopifnot(
  all(colnames(data_all) == metadata_all$samples),
  !any(duplicated(metadata_all$samples))
)

################################################################################
## SAVE
################################################################################

saveRDS(
  list(counts = data_all, metadata = metadata_all),
  file = paste0(outdir,"/proseq_combined.rds")
)
#################################################################################
data <- data_all
data_all <- NULL
metadata <- metadata_all
# correlation plot
corMT <- cor(data)
diag(corMT) <- NA

metadata <- metadata %>%
  mutate(
    Condition = factor(Condition, levels = c("ut_ut", "IFNb_1.5","IFNb_4", "IFNb_24", "IFNb_48",
                                             "IFNg_1.5", "IFNg_4", "IFNg_24", "IFNg_48")), # keep original levels
    Genotype  = factor(Genotype, levels = c("WT", "IRF2KO","IRF1KO")),
    Dataset   = factor(Dataset, levels = c("d1", "d2"))
  )

######################
design <- model.matrix(~ 0 + Condition:Genotype:Dataset, data = metadata)
pheatmap(design)
dge <- DGEList(data)
dge <- calcNormFactors(dge, method = "TMM")
keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr,]

#voom
dataVoom <- voom(dge, design=design, plot = TRUE) # insert your model matrix
# PCA
# PCA on voom-transformed data
pr_data <- file.path(outdir, "Processed_data")
if (!dir.exists(pr_data)) dir.create(pr_data, recursive = TRUE)
voom_mat <- dataVoom$E   # log2-CPM values
write_rds(voom_mat,paste0(pr_data,"/normalized_data.rds"))
# Perform PCA
pca <- prcomp(t(voom_mat), scale. = TRUE)

# Percent variance explained
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1 <- round(percentVar[1], 1)
pc2 <- round(percentVar[2], 1)

# Build PCA dataframe
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  samples = rownames(pca$x)
) %>%
  left_join(metadata, by = c("samples"))

# Relevel factors
metadata$Condition <- factor(metadata$Condition,
                             levels = c("ut_ut","IFNb_1.5","IFNb_4","IFNb_24", "IFNb_48",
                                        "IFNg_1.5", "IFNg_4", "IFNg_24", "IFNg_48"))
metadata$Genotype <- factor(metadata$Genotype, levels = c("WT", "IRF2KO","IRF1KO"))
metadata$Dataset <- factor(metadata$Dataset, levels = c("d1", "d2"))

# Color palette for Condition
my_colors <- c(
  "grey","#492050", "#82498C", "#B574C2", "#D2A9DB",   
  "#256C26","#91C392", "#4E9D4F","#C8E1C9"
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, shape = Genotype, fill = Dataset)) +
  geom_point(size = 4, stroke = 1.5) +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = c(21, 22, 24)) + # WT=circle, IRF2KO=square, IRF1KO=triangle
  scale_fill_manual(values = c("white", "black")) + # Dataset d1 = white, d2 = black
  theme_bw(base_size = 14) +
  labs(
    x = paste0("PC1 (", pc1, "%)"),
    y = paste0("PC2 (", pc2, "%)"),
    title = "PCA on normalized expression data",
    fill = "Dataset"
  )

# Save plot
qc_data <- file.path(outdir, "QC_and_basic_plots")
if (!dir.exists(qc_data)) dir.create(pr_data, recursive = TRUE)
pca <- file.path(qc_data, "PCA_on normalized_expression.pdf")
ggsave(pca)
##################################
## ============================
## CONTRAST SETUP
## ============================
####################
colnames(design) <- gsub("Dataset", "", colnames(design))
colnames(design) <- make.names(colnames(design))
conds <- levels(metadata$Condition)
conds <- conds[conds != "ut_ut"]

coef_name <- function(condition, genotype, dataset) {
  hits <- grep(
    paste0("^Condition", condition, "\\.Genotype", genotype, "\\.", dataset, "$"),
    colnames(design),
    value = TRUE
  )
  
  stopifnot(length(hits) == 1)
  hits
}
coef_name("IFNb_1.5", "IRF2KO", "d1")
coef_name("ut_ut",    "WT",     "d2")

conds <- levels(metadata$Condition)
conds <- conds[conds != "ut_ut"]


## ============================
## WITHIN-DATASET INTERACTIONS
## ============================
baseline_IRF1KO_vs_WT <- paste0(
  coef_name("ut_ut", "IRF1KO", "d2"),
  " - ",
  coef_name("ut_ut", "WT", "d2")
)
names(baseline_IRF1KO_vs_WT) <- "IRF1KO_ut_vs_WT"
baseline_IRF1KO_vs_WT
baseline_IRF2KO_vs_WT <- paste0(
  coef_name("ut_ut", "IRF2KO", "d1"),
  " - ",
  coef_name("ut_ut", "WT", "d1")
)

names(baseline_IRF2KO_vs_WT) <- "IRF2KO_ut_vs_WT"

# d1: IRF2KO vs WT
interaction_d1 <- sapply(conds, function(cond) {
  paste0(
    "(", coef_name(cond, "IRF2KO", "d1"), " - ", coef_name("ut_ut", "IRF2KO", "d1"), ") - ",
    "(", coef_name(cond, "WT",     "d1"), " - ", coef_name("ut_ut", "WT",     "d1"), ")"
  )
})
names(interaction_d1) <- paste0("Interaction_d1_IRF2KO_vs_WT_", conds)

# d2: IRF1KO vs WT
interaction_d2 <- sapply(conds, function(cond) {
  paste0(
    "(", coef_name(cond, "IRF1KO", "d2"), " - ", coef_name("ut_ut", "IRF1KO", "d2"), ") - ",
    "(", coef_name(cond, "WT",     "d2"), " - ", coef_name("ut_ut", "WT",     "d2"), ")"
  )
})
names(interaction_d2) <- paste0("Interaction_d2_IRF1KO_vs_WT_", conds)

## ============================
## CROSS-DATASET INTERACTIONS
## ============================

cross_dataset <- sapply(conds, function(cond) {
  paste0(
    "(", coef_name(cond, "IRF2KO", "d1"), " - ", coef_name("ut_ut", "IRF2KO", "d1"), ") - ",
    "(", coef_name(cond, "WT",     "d1"), " - ", coef_name("ut_ut", "WT",     "d1"), ") - ",
    "(", coef_name(cond, "IRF1KO", "d2"), " - ", coef_name("ut_ut", "IRF1KO", "d2"), ") + ",
    "(", coef_name(cond, "WT",     "d2"), " - ", coef_name("ut_ut", "WT",     "d2"), ")"
  )
})
names(cross_dataset) <- paste0("IRF2KO_vs_IRF1KO_", conds)
ut_KO <- paste0(
"(",coef_name("ut_ut",    "IRF2KO",     "d1"), " - ", coef_name("ut_ut", "WT", "d1"), ") - ",
"(",coef_name("ut_ut",    "IRF1KO",     "d2"), " - ", coef_name("ut_ut", "WT", "d2"), ")"
)
names(ut_KO) <-   "IRF2KO_vs_IRF1KO_ut"

## ============================
## COMBINE ALL CONTRASTS
## ============================

all_contrasts <- c(interaction_d1, interaction_d2, cross_dataset, ut_KO,
                   baseline_IRF2KO_vs_WT,baseline_IRF1KO_vs_WT)

contrast_matrix <- makeContrasts(
  contrasts = all_contrasts,
  levels = design
)
colnames(contrast_matrix) <- names(all_contrasts)
## ============================
## LIMMA PIPELINE
## ============================

dge <- DGEList(data)
dge <- calcNormFactors(dge, method = "TMM")
keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr,]

#voom
dataVoom <- voom(dge, design=design, plot = TRUE)
fit  <- lmFit(dataVoom, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

## ============================
## COLLECT RESULTS
## ============================
contrast_names <- colnames(contrast_matrix)
pdf(paste0(outdir,"Contrast.pdf"), h=15, w = 12)
pheatmap(contrast_matrix)
dev.off()
limmaRes <- lapply(colnames(contrast_matrix), function(cn) {
  topTable(fit2, coef = cn, number = Inf) %>%
    rownames_to_column("Gene") %>%
    mutate(Contrast = cn)%>%
    mutate(group = case_when(
             logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
             logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
             TRUE ~ "n.s"
           ))
}) %>%
  bind_rows()
t1 <- limmaRes %>%
  count(Contrast, group)

##################
#heatmap--------------
##################
#selection here-
#significant interaction

unique(limmaRes$Contrast)
grep("^IRF", unique(limmaRes$Contrast), value = T)
## ============================
## EXTRACT OPPOSITE-TREND GENES--------------NOWT INDUCTION CONSIDERED
## ============================

## 1. Split within-dataset interaction results

# res_d1 <- limmaRes %>%
#   filter(Contrast %in% c("IRF2KO_ut_vs_WT",grep("^Interaction_d1_IRF2KO_vs_WT_",
#                                                 Contrast, value = T))) %>%
#   mutate(
#     Condition = sub("^Interaction_d1_IRF2KO_vs_WT_", "", Contrast)
#   )%>%
#   mutate(Condition = gsub("IRF2KO_ut_vs_WT","ut", Condition))
# 
# res_d2 <- limmaRes %>%
#   filter(Contrast %in% c("IRF1KO_ut_vs_WT",grep("^Interaction_d2_IRF1KO_vs_WT_",
#                                                 Contrast,value = T))) %>%
#   mutate(
#     Condition = sub("^Interaction_d2_IRF1KO_vs_WT_", "", Contrast)
#   )%>%
#   mutate(Condition = gsub("IRF1KO_ut_vs_WT","ut", Condition))
# 
# ## 2. Join IRF2KO (d1) and IRF1KO (d2) effects per gene + condition
# res_joined <- res_d1 %>%
#   dplyr::select(
#     Gene,
#     Condition,
#     logFC_IRF2KO = logFC,
#     P_d1     = P.Value,
#     adjP_IRF2KO  = adj.P.Val
#   ) %>%
#   inner_join(
#     res_d2 %>%
#       dplyr::select(
#         Gene,
#         Condition,
#         logFC_IRF1KO = logFC,
#         P_d2     = P.Value,
#         adjP_IRF1KO  = adj.P.Val
#       ),
#     by = c("Gene", "Condition")
#   )%>%
#   mutate(group = ifelse(abs(logFC_IRF2KO) > 0.5 & 
#                           abs(logFC_IRF1KO) > 0.5 &
#                           adjP_IRF2KO < 0.05 &
#                           adjP_IRF1KO < 0.05, "both",
#                         ifelse(abs(logFC_IRF1KO) > 0.5 & adjP_IRF1KO < 0.05,"IRF1", 
#                                ifelse(abs(logFC_IRF2KO) > 0.5 & adjP_IRF2KO < 0.05,"IRF2" ,"ns"))))
# 
# unique(res_joined$Condition)
# res_joined$Condition <- factor(res_joined$Condition, levels = c("ut", "IFNb_1.5","IFNb_4", "IFNb_24", "IFNb_48",
#                                                      "IFNg_1.5", "IFNg_4", "IFNg_24", "IFNg_48"))
# 
# ggplot(res_joined, aes(x = logFC_IRF2KO, y = logFC_IRF1KO)) +
#   
#   # 1️⃣ Background hex for ns
#   geom_hex(
#     data = dplyr::filter(res_joined, group == "ns"),
#     bins = 50,
#     alpha = 0.6
#   ) +
#   scale_fill_gradient(low = "grey", high = "grey40", guide = "none") +
#   
#   # 2️⃣ IRF1 & IRF2 first
#   geom_point(
#     data = dplyr::filter(res_joined, group %in% c("IRF1", "IRF2")),
#     aes(color = group),
#     alpha = 0.7,
#     size = 1.2
#   ) +
#   
#   # 3️⃣ BOTH last (strictly on top)
#   geom_point(
#     data = dplyr::filter(res_joined, group == "both"),
#     aes(color = group),
#     alpha = 0.9,
#     size = 1.8
#   ) +
#   
#   scale_color_manual(
#     values = c(
#       IRF2 = "#1C3C13",
#       IRF1 = "#4A3563",
#       both = "#A24C06"
#     ),
#     name = "Significant in"
#   ) +
#   
#   coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) +
#   facet_wrap(~Condition) +
#   optimized_theme_fig() +
#   theme(panel.grid.major = NULL)
# #save
# volcano_dir <- file.path(outdir, "Comparrison_volcano_plots")
# if (!dir.exists(volcano_dir)) dir.create(volcano_dir, recursive = TRUE)
# ggsave(paste0(volcano_dir,"/significant_genes_IRF1_IRF2_both.pdf"))
#################################
#Considering WT induction----------------

WT <- read_rds(here("Results/IRF2_01_DE/Processed_data/","/limmaRes.rds"))
unique(WT$coef)
WT <- WT |>
  filter(!(grepl("IRF2",WT$coef)))
head(WT)

WT <- WT |>
  mutate(Condition = coef)|>
  dplyr::select(
    Gene = gene.name,
    Condition,
    logFC_WT = logFC,
    P_WT     = P.Value,
    adjP_WT  = adj.P.Val
  )
## 1. Split within-dataset interaction results

res_d1 <- limmaRes %>%
  filter(Contrast %in% c("IRF2KO_ut_vs_WT",grep("^Interaction_d1_IRF2KO_vs_WT_",
                                                Contrast, value = T))) %>%
  mutate(
    Condition = sub("^Interaction_d1_IRF2KO_vs_WT_", "", Contrast)
  )%>%
  mutate(Condition = gsub("IRF2KO_ut_vs_WT","ut", Condition))

res_d2 <- limmaRes %>%
  filter(Contrast %in% c("IRF1KO_ut_vs_WT",grep("^Interaction_d2_IRF1KO_vs_WT_",
                                                Contrast,value = T))) %>%
  mutate(
    Condition = sub("^Interaction_d2_IRF1KO_vs_WT_", "", Contrast)
  )%>%
  mutate(Condition = gsub("IRF1KO_ut_vs_WT","ut", Condition))

## 2. Join IRF2KO (d1) and IRF1KO (d2) effects per gene + condition
res_joined_MOD <- res_d1 %>%
  dplyr::select(
    Gene,
    Condition,
    logFC_IRF2KO = logFC,
    P_d1     = P.Value,
    adjP_IRF2KO  = adj.P.Val
  ) %>%
  inner_join(
    res_d2 %>%
      dplyr::select(
        Gene,
        Condition,
        logFC_IRF1KO = logFC,
        P_d2     = P.Value,
        adjP_IRF1KO  = adj.P.Val
      ),
    by = c("Gene", "Condition")
  )%>%
  inner_join(
    WT %>%
      dplyr::select(
        Gene,
        Condition,
        logFC_WT,
        P_WT,
        adjP_WT
      ),
    by = c("Gene", "Condition"))|>
  mutate(group = ifelse(abs(logFC_IRF2KO) > 0.5 & 
                          abs(logFC_IRF1KO) > 0.5 &
                          abs(logFC_WT) >0.5 &
                          adjP_WT < 0.05 &
                          adjP_IRF2KO < 0.05 &
                          adjP_IRF1KO < 0.05, "both",
                        ifelse(abs(logFC_IRF1KO) > 0.5 & adjP_IRF1KO < 0.05 & adjP_WT < 0.05 & abs(logFC_WT)> 0.5,"IRF1", 
                               ifelse(abs(logFC_IRF2KO) > 0.5 & adjP_IRF2KO < 0.05 & adjP_WT < 0.05 & abs(logFC_WT)> 0.5,"IRF2" ,"ns"))))

unique(res_joined_MOD$Condition)
res_joined_MOD$Condition <- factor(res_joined_MOD$Condition, levels = c("ut", "IFNb_1.5","IFNb_4", "IFNb_24", "IFNb_48",
                                                                "IFNg_1.5", "IFNg_4", "IFNg_24", "IFNg_48"))


write_rds(res_joined_MOD,paste0(pr_data,"/res_joined_IRF1_IRF2_significant_genes.rds"))
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
res_both <- res_joined_MOD %>% filter(group != "ns")%>%
  filter(group == "both")

#function to plot--------------
plot_condition_heatmap <- function(cond_name,
                                   res_both,
                                   expr_matrix,
                                   metadata,
                                   outdir,
                                   label) {
  
  message("Processing: ", cond_name)
  
  # ---- 1. Get genes significant in this condition ----
  genes_cond <- res_both %>%
    filter(Condition == cond_name) %>%
    pull(Gene) %>%
    unique()
  
  if(length(genes_cond) == 0) {
    message("No genes for ", cond_name)
    return(NULL)
  }
  
  # ---- 2. Get UT + this condition samples ----
  metadata_sub <- metadata %>%
    filter(Condition %in% c("ut_ut", cond_name))
  
  # Order samples: Condition → Dataset → Genotype → Replicate
  metadata_sub <- metadata_sub %>%
    mutate(
      Condition = factor(Condition, levels = c("ut_ut", cond_name)),
      Dataset   = factor(Dataset, levels = c("d1","d2")),
      Genotype  = factor(Genotype, levels = c("WT","IRF2KO","IRF1KO"))
    ) %>%
    arrange(
      Dataset,
      Genotype,
      Condition,
      Replicate
    )
  
  
  samples_use <- metadata_sub$samples
  
  # ---- 3. Subset expression matrix ----
  mat <- expr_matrix[genes_cond, samples_use, drop = FALSE]
  
  # Remove genes not found
  mat <- mat[complete.cases(mat), , drop = FALSE]
  
  if(nrow(mat) == 0) {
    message("Matrix empty for ", cond_name)
    return(NULL)
  }
  mat_z <- mat
  
  
  # ---- Define condition colors dynamically ----
  cond_colors <- c("ut_ut" = "#ABC5BB")
  cond_colors[cond_name] <- "#3C549A"
  col_ha <- HeatmapAnnotation(
    Condition = metadata_sub$Condition,
    Dataset   = metadata_sub$Dataset,
    Genotype  = metadata_sub$Genotype,
    col = list(
      Condition = cond_colors,
      Genotype = c(
        WT      = "grey60",
        IRF2KO  = "#1C3C13",
        IRF1KO  = "#4A3563"
      ),
      Dataset = c(
        d1 = "#6EC187",
        d2 = "#966AC3"
      )
    )
  )
  
  
  ht <- Heatmap(
    mat_z,
    name = "label",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = T,
    show_row_dend = FALSE,
    top_annotation = col_ha,
    column_title = paste("Genes significant in BOTH —", cond_name),
    col = colorRamp2(c(-2, 0, 2),
                     c("#4C889C", "white", "#D0154E")),
    column_split = metadata_sub$Dataset,   # <-- split by dataset
    column_gap = unit(4, "mm")
  )
  
  # ---- 6. Dynamic height ----
  pdf_height <- max(5, min(40, 0.25 * nrow(mat_z)))
  
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  pdf(file.path(outdir,
                paste0("Heatmap_both_", cond_name, ".pdf")),
      ,
      width = 8,
      height = pdf_height)
  
  draw(ht)
  
  dev.off()
}
cond_to_plot <- NULL
conditions_to_plot <-NULL
cond_to_plot <- unique(res_both$Condition)
conditions_to_plot <- c(
  "IFNb_1.5", "IFNb_4",   "IFNg_4",
  "IFNg_24",  "IFNg_48"
)

#expr with zscore---------
#########original_expr_zscore------------
expr_matrix <- dataVoom$E

IRF1 <- expr_matrix[, rownames(metadata2)]
IRF2 <- expr_matrix[, rownames(metadata1)]

IRF1 <- t(scale(t(IRF1)))
IRF2 <- t(scale(t(IRF2))) 
expr_matrix <- cbind(IRF1,IRF2)

for(cond in conditions_to_plot){
  plot_condition_heatmap(
    cond_name = cond,
    res_both = res_both,
    expr_matrix = expr_matrix,   # ✅ correct
    metadata = metadata,
    outdir = "IRF1_vs_IRF2_combinedTMM/Expr_zscore"
  )
}



#remove batch effect--------------
bio_design <- model.matrix(
  ~ 0 + Condition:Genotype,
  data = metadata
)

colnames(bio_design) <- make.names(colnames(bio_design))

#batch effect removal step
expr_matrix <- removeBatchEffect(
  dataVoom$E,
  batch  = metadata$Dataset,
  design = bio_design
)
#
IRF1 <- expr_matrix[, rownames(metadata2)]
IRF2 <- expr_matrix[, rownames(metadata1)]

IRF1 <- t(scale(t(IRF1)))
IRF2 <- t(scale(t(IRF2))) 
expr_matrix <- cbind(IRF1,IRF2)
for(cond in conditions_to_plot){
  plot_condition_heatmap(
    cond_name = cond,
    res_both = res_both,
    expr_matrix = expr_matrix,   # ✅ correct
    metadata = metadata,
    outdir = paste0(outdir,"/Expression_plots_per_condition",
                    "/remove_batch_expr_zscore")
  )
}




# ###########
# #original expr but divide by mean instead of sd--------
# expr_matrix <- dataVoom$E
# mean_scale <- function(mat) {
#   t(apply(mat, 1, function(x) {
#     m <- mean(x, na.rm = TRUE)
#     if (m == 0) return(rep(0, length(x)))
#     (x - m) / m
#   }))
# }
# IRF1 <- expr_matrix[, rownames(metadata2)]
# IRF2 <- expr_matrix[, rownames(metadata1)]
# IRF1 <- mean_scale(IRF1)
# IRF2 <- mean_scale(IRF2) 
# expr_matrix <- cbind(IRF1,IRF2)
# 
# 
# for(cond in conditions_to_plot){
#   plot_condition_heatmap(
#     cond_name = cond,
#     res_both = res_both,
#     expr_matrix = expr_matrix,   # ✅ correct
#     metadata = metadata,
#     outdir = "IRF1_vs_IRF2_combinedTMM/Expr_mean_division",
#     label = "rel_expr"
#   )
# }
# ######################
# #mean---on batch removed---------
# bio_design <- model.matrix(
#   ~ 0 + Condition:Genotype,
#   data = metadata
# )
# 
# colnames(bio_design) <- make.names(colnames(bio_design))
# 
# #batch effect removal step
# expr_matrix <- removeBatchEffect(
#   dataVoom$E,
#   batch  = metadata$Dataset,
#   design = bio_design
# )
# 
# mean_scale <- function(mat) {
#   t(apply(mat, 1, function(x) {
#     m <- mean(x, na.rm = TRUE)
#     if (m == 0) return(rep(0, length(x)))
#     (x - m) / m
#   }))
# }
# IRF1 <- expr_matrix[, rownames(metadata2)]
# IRF2 <- expr_matrix[, rownames(metadata1)]
# IRF1 <- mean_scale(IRF1)
# IRF2 <- mean_scale(IRF2) 
# expr_matrix <- cbind(IRF1,IRF2)
# 
# 
# for(cond in conditions_to_plot){
#   plot_condition_heatmap(
#     cond_name = cond,
#     res_both = res_both,
#     expr_matrix = expr_matrix,   # ✅ correct
#     metadata = metadata,
#     outdir = "IRF1_vs_IRF2_combinedTMM/batch_remove_expr_mean_division",
#     label = "rel_expr"
#   )
# }
# ###################
# 
# expr_matrix <- removeBatchEffect(
#   dataVoom$E,
#   batch  = metadata$Dataset,
#   design = bio_design
# )
# expr_matrix <- t(scale(t(expr_matrix))) 
# 
# for(cond in conditions_to_plot){
#   plot_condition_heatmap(
#     cond_name = cond,
#     res_both = res_both,
#     expr_matrix = expr_matrix,   # ✅ correct
#     metadata = metadata,
#     outdir = "IRF1_vs_IRF2_combinedTMM/zscore_across_dataset",
#     label = "rel_expr"
#   )
# }
# 
