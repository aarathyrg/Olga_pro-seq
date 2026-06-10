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
library(biomaRt)
library(dplyr)
source("Ag_optimized_theme.R")
################################################################################
outdir <- here("Results/IRF2_04_reanalyze_IRF1/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

pr_data <- file.path(outdir, "Processed_data")
if (!dir.exists(pr_data)) dir.create(pr_data, recursive = TRUE)
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


## ---- Align genes ----
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


data2 <- data2[common_genes, ]

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


rownames(metadata2) <- metadata2$samples

if (!dir.exists(pr_data)) dir.create(pr_data, recursive = TRUE)
outfile <- file.path(pr_data, "metadata.rds")
metadata2 %>% write_rds(outfile)
################################################################################
## QC: SAMPLE CORRELATION
################################################################################
corMT <- cor(data2)
diag(corMT) <- NA

qcdir <- file.path(outdir, "QC_and_basic_plots")
if (!dir.exists(qcdir)) dir.create(qcdir, recursive = TRUE)

pdf(file.path(qcdir, "sample_correlation_IRF1.pdf"), width = 25, height = 25)
Heatmap(corMT,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12))
dev.off()

################################################################################
## DESIGN MATRIX (IDENTICAL STRUCTURE TO IRF2)
################################################################################
design <- model.matrix(~Condition * Genotype, data = metadata2)

dge <- DGEList(data2)
dge <- calcNormFactors(dge, method = "TMM")

keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr, ]

################################################################################
## VOOM TRANSFORMATION
################################################################################
dataVoom <- voom(dge, design, plot = TRUE)
outfile <- file.path(pr_data, "dataVoom.rds")
dataVoom %>% write_rds(outfile)


################################################################################
## PCA
################################################################################
voom_mat <- dataVoom$E

pca <- prcomp(t(voom_mat), scale. = TRUE)

percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  samples = rownames(pca$x)
) %>%
  left_join(metadata2, by = c("samples"))
metadata2$samples
ggplot(pca_df, aes(PC1, PC2, color = Condition, shape = Genotype)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "IRF1 PCA",
       x = paste0("PC1 (", round(percentVar[1],1), "%)"),
       y = paste0("PC2 (", round(percentVar[2],1), "%)"))

ggsave(file.path(qcdir, "PCA_IRF1.pdf"))

################################################################################
## LIMMA MODEL
################################################################################
limmaFit <- lmFit(dataVoom, design)
limmaFit <- eBayes(limmaFit)

################################################################################
## EXTRACT COEFFICIENTS
################################################################################
limmaRes <- map_dfr(colnames(coef(limmaFit)), function(coefx) {
  topTable(limmaFit, coef = coefx, number = Inf) %>%
    rownames_to_column("genes") %>%
    filter(coefx != "(Intercept)") %>%
    mutate(
      coef = coefx,
      group = case_when(
        logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
        logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
        TRUE ~ "n.s"
      )
    )
})

################################################################################
## CONTRASTS (IRF1 VERSION - SAME STRUCTURE AS IRF2)
################################################################################
treatments <- c("IFNb_1.5","IFNb_4","IFNb_24","IFNb_48",
                "IFNg_1.5","IFNg_4","IFNg_24","IFNg_48")

colnames(limmaFit$coefficients) <- gsub(":", ".", colnames(limmaFit$coefficients))

contrast_list <- sapply(treatments, function(trt) {
  main <- paste0("Condition", trt)
  inter <- paste0("Condition", trt, ".GenotypeIRF1KO")
  paste0(main, "+", inter)
})

contrast_matrix <- makeContrasts(
  contrasts = contrast_list,
  levels = colnames(limmaFit$coefficients)
)

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
limmaRes <- bind_rows(limmaRes, limmaRes_contrasts)

limmaRes <- limmaRes %>%
  mutate(
    # Remove prefixes
    coef = str_replace_all(coef, "Condition|Genotype", ""),
    coef = str_replace(coef, "IRF1-/-", "IRF1KO"),
    coef = if_else(
      str_detect(coef, "\\+"),
      str_replace(coef, "\\+.*", "_IRF1KO"), coef ),
    # Handle interactions (contains ":")
    coef = if_else(
      str_detect(coef, ":"),
      # Replace b/g at start with IFNb/IFNg, replace : with _, append _Interaction
      coef %>% str_replace(":", "_") %>%
        paste0("_Interaction"),
      coef
    ))
unique(limmaRes$coef)
################################################################################
## GENE ANNOTATION
################################################################################
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = unique(limmaRes$genes),
  mart = ensembl
)

limmaRes$gene.name <- gene_info$mgi_symbol[
  match(limmaRes$genes, gene_info$ensembl_gene_id)
]


outfile <- file.path(pr_data, "limmaRes_IRF1.rds")
limmaRes %>% write_rds(outfile)
# Split limmaRes into a list of tibbles, one per coefficient
limmaRes_list <- limmaRes %>%
  group_by(coef) %>%
  group_split() %>%
  setNames(unique(limmaRes$coef))

# Create folder if it doesn't exist
# Split by coefficient (CORRECT & SAFE)
limmaRes_list <- split(limmaRes, limmaRes$coef)

# Create output directory if needed
res_dir <- file.path(pr_data, "Table_of_results_IRF1")


if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}

# Save limma results per coefficient
for (coef_name in names(limmaRes_list)) {
  write.csv(
    limmaRes_list[[coef_name]],
    file = file.path(
      res_dir,
      paste0("limma_results_IRF1", coef_name, ".csv")
    ),
    row.names = FALSE
  )
}



#Significant results----------------------
limmaRessig <- limmaRes %>%
  filter(adj.P.Val < 0.05, abs(logFC)> 1)
# Split SIGNIFICANT results by coefficient
limmaResSig_list <- split(limmaRessig, limmaRessig$coef)

# Create output directory if needed
sig_res_dir <- file.path(pr_data, "Table_of_SIGNIFICANT_results_IRF1")


if (!dir.exists(sig_res_dir)) {
  dir.create(sig_res_dir)
}

# Save significant limma results per coefficient
for (coef_name in names(limmaResSig_list)) {
  write.csv(
    limmaResSig_list[[coef_name]],
    file = file.path(
      sig_res_dir,
      paste0("limma_SIGNIFICANT_results_IRF1", coef_name, ".csv")
    ),
    row.names = FALSE
  )
}
#############################################

