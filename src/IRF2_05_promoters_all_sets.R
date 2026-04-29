library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)

# Parameters
lfc_cutoff  <- 0.5
padj_cutoff <- 0.05
upstream    <- 2000
downstream  <- 200
#####
outdir <- here("Results/IRF2_05_promoters_all_sets/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
#Input
res_joined <- read_rds(here("Results/IRF2_04_comparison_with_IRF1_dataset/Processed_data/",
                            "res_joined_IRF1_IRF2_significant_genes.rds"))


# Function: IRF2-repressed, IRF1-activated genes
get_IRF2_rep_IRF1_act <- function(res_joined, condition,
                                  lfc_cutoff = 0.5,
                                  padj_cutoff = 0.05) {
  res_joined %>%
    filter(Condition == condition) %>%
    filter(
      logFC_IRF2KO >= lfc_cutoff,   # up in IRF2KO → repressed by IRF2
      logFC_IRF1KO <= -lfc_cutoff,  # down in IRF1KO → activated by IRF1
      adjP_IRF2KO < padj_cutoff,
      adjP_IRF1KO < padj_cutoff
    ) %>%
    pull(Gene) %>%
    unique()
}


# Define background sets
get_IRF2_act <- function(res_joined, condition, lfc_cutoff = 0.5, padj_cutoff = 0.05) {
  res_joined %>%
    filter(Condition == condition) %>%
    filter(logFC_IRF2KO <= -lfc_cutoff, adjP_IRF2KO < padj_cutoff) %>%
    pull(Gene) %>%
    unique()
}

get_not_IRF2_reg <- function(res_joined, condition) {
  res_joined %>%
    filter(Condition == condition) %>%
    filter(adjP_IRF2KO > 0.2) %>%
    pull(Gene) %>%
    unique()
}

# Conditions
conditions_to_plot <- c(
  "IFNb_1.5", "IFNb_4", 
  "IFNg_4", "IFNg_24", "IFNg_48"
)

# Map gene symbols → Entrez IDs
all_genes <- unique(c(
  unlist(lapply(conditions_to_plot, function(cond) get_IRF2_rep_IRF1_act(res_joined, cond))),
  unlist(lapply(conditions_to_plot, function(cond) get_IRF2_act(res_joined, cond))),
  unlist(lapply(conditions_to_plot, function(cond) get_not_IRF2_reg(res_joined, cond)))
))

gene_map <- mapIds(
  org.Mm.eg.db,
  keys = all_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Function: extract promoters and save BED + FASTA
extract_promoters <- function(gene_symbols, outdir, label,
                              upstream = 2000, downstream = 200) {
  entrez_ids <- na.omit(gene_map[gene_symbols])
  txdb_genes <- names(genes(txdb))
  entrez_ids <- entrez_ids[entrez_ids %in% txdb_genes]
  
  if(length(entrez_ids) == 0){
    message("No genes to process for ", label)
    return(NULL)
  }
  
  # Get promoter GRanges
  grl <- genes(txdb, single.strand.genes.only = FALSE)[entrez_ids]
  prom_grl <- promoters(grl, upstream = upstream, downstream = downstream)
  prom_gr <- unlist(prom_grl)
  
  # Export BED
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  bed_file <- file.path(outdir, paste0(label, ".bed"))
  rtracklayer::export(prom_gr, bed_file)
  
  # Export FASTA
  fasta_file <- file.path(outdir, paste0(label, ".fa"))
  seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, prom_gr)
  names(seqs) <- prom_gr$gene_id
  Biostrings::writeXStringSet(seqs, fasta_file)
  
  message("Saved ", length(prom_gr), " promoters for ", label)
  return(prom_grl)
}

# Process foreground + backgrounds


for(cond in conditions_to_plot){
  message("Processing condition: ", cond)
  
  fg_genes <- get_IRF2_rep_IRF1_act(res_joined, cond)
  bg1_genes <- get_IRF2_act(res_joined, cond)
  bg2_genes <- get_not_IRF2_reg(res_joined, cond)
  
  extract_promoters(fg_genes, outdir, paste0(cond, "_IRF2rep_IRF1act"))
  extract_promoters(bg1_genes, outdir, paste0(cond, "_IRF2act_bg"))
  extract_promoters(bg2_genes, outdir, paste0(cond, "_not_IRF2_bg"))
}
