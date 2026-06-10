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
outdir <- here("Results/IRF2_06_IRF2_promoters_all_sets/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Input
limmaRes <- read_rds(here("Results/IRF2_01_DE/Processed_data/limmaRes.rds"))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

############################
# GENE MAP (your style, fixed once)
############################

all_genes <- unique(limmaRes$gene.name)

gene_map <- mapIds(
  org.Mm.eg.db,
  keys = all_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

############################
# gene sets
############################

IRF2_repressed_ut <- limmaRes |>
  filter(coef == "IRF2KO") |>
  filter(adj.P.Val < 0.05 & logFC > 1) |>
  pull(gene.name) |>
  unique()

IRF2_activated_ut <- limmaRes |>
  filter(coef == "IRF2KO") |>
  filter(adj.P.Val < 0.05 & logFC < -1) |>
  pull(gene.name) |>
  unique()

IRF2_independent <- limmaRes |>
  filter(coef == "IRF2KO") |>
  filter(adj.P.Val > 0.1) |>
  pull(gene.name) |>
  unique()

############################
# FUNCTION (fixed gmap usage + TxDb indexing)
############################

extract_promoters <- function(gene_symbols, outdir, label,
                              upstream = 2000, downstream = 200) {
  
  # use GLOBAL gene_map (your style)
  # convert SYMBOL → ENTREZ safely
  entrez_ids <- gene_map[gene_symbols]
  entrez_ids <- na.omit(entrez_ids)
  entrez_ids <- unique(entrez_ids)
  
  txdb_genes <- names(genes(txdb))
  entrez_ids <- entrez_ids[entrez_ids %in% txdb_genes]
  
  if (length(entrez_ids) == 0) {
    message("No genes to process for ", label)
    return(NULL)
  }
  
  # FIX: proper TxDb subsetting (no GRL list indexing)
  tx_genes <- genes(txdb)
  tx_genes <- tx_genes[names(tx_genes) %in% entrez_ids]
  
  # promoters
  prom_gr <- promoters(tx_genes,
                       upstream = upstream,
                       downstream = downstream)
  
  # Export BED
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  bed_file <- file.path(outdir, paste0(label, ".bed"))
  rtracklayer::export(prom_gr, bed_file)
  
  # Export FASTA
  fasta_file <- file.path(outdir, paste0(label, ".fa"))
  seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, prom_gr)
  names(seqs) <- names(prom_gr)
  
  Biostrings::writeXStringSet(seqs, fasta_file)
  
  message("Saved ", length(prom_gr), " promoters for ", label)
  
  return(prom_gr)
}

############################
############################
# PROCESS PROMOTERS
############################

# IRF2 repressed (foreground)
prom_repressed <- extract_promoters(
  IRF2_repressed_ut,
  outdir,
  "IRF2_repressed_ut"
)

# IRF2 activated (foreground)
prom_activated <- extract_promoters(
  IRF2_activated_ut,
  outdir,
  "IRF2_activated_ut"
)

# non-regulated background
prom_background <- extract_promoters(
  IRF2_independent,
  outdir,
  "IRF2_independent_bg"
)
