coefx <- "IRF2KO"

### ---- 1. Get significant genes ----
sig_genes <- limmaRessig %>%
  filter(coef == coefx) %>%
  pull(genes) %>% unique()
if (length(sig_genes) == 0) stop("No sig genes")

### ---- 2. Samples ----
relevant_samples <- get_relevant_samples(coefx, metadata)

### ---- 3. Expression matrix ----
expr_df <- scale(dataVoom$E[sig_genes, relevant_samples]) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

expr_long <- expr_df %>%
  pivot_longer(cols = -gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(metadata |> filter(Sample %in% relevant_samples),
            by = c("Sample")) %>%
  mutate(gene_name = limmaRessig$gene.name[match(gene, limmaRessig$genes)]) %>%
  drop_na(gene_name)

### ---- 4. Sample ordering ----
sample_order <- {
  if (coefx == "IRF2KO") {
    c(
      metadata |> filter(Genotype=="WT",     Time=="ut") |> pull(samples),
      metadata |> filter(Genotype=="IRF2KO", Time=="ut") |> pull(samples)
    )
  } else {
    prefix <- sub("_IRF2KO_Interaction", "", coefx)
    c(
      metadata |> filter(Genotype=="WT",     Time=="ut") |> pull(samples),
      metadata |> filter(Genotype=="WT",     grepl(paste0("^", prefix), Condition)) |> pull(samples),
      metadata |> filter(Genotype=="IRF2KO", Time=="ut") |> pull(samples),
      metadata |> filter(Genotype=="IRF2KO", grepl(paste0("^", prefix), Condition)) |> pull(samples)
    ) |> unique()
  }
}

expr_long$Sample <- factor(expr_long$Sample, levels = sample_order)

### ---- 5. Heatmap matrix ----
heat_mat <- expr_long %>%
  select(gene_name, Sample, Expression) %>%
  pivot_wider(names_from = Sample, values_from = Expression) %>%
  tibble::column_to_rownames("gene_name") %>%
  as.matrix()

heat_mat <- heat_mat[, sample_order[sample_order %in% colnames(heat_mat)], drop = FALSE]

### ---- color function ----
col_fun <- circlize::colorRamp2(
  c(min(heat_mat, na.rm=TRUE), 0, max(heat_mat, na.rm=TRUE)),
  c("#4C889C","white","#D0154E")
)

### ---- 6. Heatmap ----
ht <- Heatmap(
  heat_mat,
  name = "Expression",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  show_row_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_title = paste("Significant genes:", coefx),
  na_col = "grey"
)

ht_drawn <- draw(ht)
ht_grob  <- grid.grabExpr(draw(ht))

ordered_genes <- rownames(heat_mat)[row_order(ht_drawn)]
n_genes <- length(ordered_genes)
pdf_height <- max(6, min(40, 6 + 0.25 * n_genes))

### ---- 7. Dot plot ----
needed <- get_related_coefs(coefx)

dot_df <- limmaRes %>%
  filter(genes %in% sig_genes, coef %in% needed) %>%
  mutate(
    coef      = factor(coef, levels = needed),
    gene.name = factor(gene.name, levels = rev(ordered_genes))
  ) %>%
  drop_na()

p <- ggplot(dot_df,
            aes(x = coef, y = gene.name,
                size = pmin(3, -log10(adj.P.Val)),
                fill = pmax(pmin(logFC, 2), -2))) +
  geom_point(shape=21, color="black") +
  scale_size_continuous(range = c(0, 3),
                        name = TeX("$-\\log_{10}(p_{adj})$")) +
  scale_fill_gradient2(high="#D0154E", mid="white", low="#4C889C") +
  labs(title = paste("Dot plot for", coefx)) +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
        axis.text.y = element_text(size=7),
        axis.ticks.y = element_blank())

### ---- 8. Save outputs ----
dir.create("Expression_plots_and_logFC_per_comparison", showWarnings = FALSE)

file_combined <- file.path("Expression_plots_and_logFC_per_comparison",
                           paste0(coefx, "_heatmap_dotplot.pdf"))

combined <- wrap_elements(full = ht_grob) | p +
  plot_layout(widths = c(1,1.8))

ggsave(file_combined, combined, width = 10, height = pdf_height, limitsize = FALSE)
