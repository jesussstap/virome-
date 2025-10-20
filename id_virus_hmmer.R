library(dplyr)
library(tidyr)
library(pheatmap)
library(tools)

# --- 1. Leer archivos HMMER ---
ruta_hmmer <- "/home/jesus/CH_22/metagenoma/resultados_hmmer/circulares/Prokka_on_data_12:_faa.fasta.tabular"
archivos <- list.files(ruta_hmmer, pattern = "\\.tabular$", full.names = TRUE)

leer_hmmer <- function(f) {
  nlineas <- sum(!grepl("^#", readLines(f)))
  if (nlineas == 0) return(NULL)
  
  df <- read.table(f, comment.char = "#", fill = TRUE, stringsAsFactors = FALSE)
  colnames(df) <- c("target_name", "accession", "tlen", 
                    "query_name", "query_accession", "qlen",
                    "full_seq_Evalue", "full_seq_score", "full_seq_bias",
                    "dom_num", "dom_of", "c_Evalue", "i_Evalue",
                    "dom_score", "dom_bias",
                    "hmm_from", "hmm_to", "ali_from", "ali_to", 
                    "env_from", "env_to", "acc", "description")
  df$sample <- file_path_sans_ext(basename(f))
  return(df)
}

lista_res <- lapply(archivos, leer_hmmer)
hmmer_res <- do.call(rbind, lista_res[!sapply(lista_res, is.null)])

# --- 2. Leer tabla de anotaciones ---
vog_annot <- read.delim("/home/jesus/CH_22/metagenoma/vog_hmm/vog.lca.tsv",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(vog_annot)[1] <- "target_name"

# --- 3. Unir resultados de HMMER con anotaciones ---
hmmer_annot <- hmmer_res %>%
  left_join(vog_annot, by = "target_name")

x# --- 4. Preparar tabla para heatmap usando taxonomía ---
tax_col <- "LastCommonAncestor_Name"  # categoría para el heatmap

# Extraer solo el nivel más bajo de LastCommonAncestor_Name
hmmer_annot <- hmmer_annot %>%
  filter(!is.na(LastCommonAncestor_Name)) %>%
  mutate(TaxonLevel = sapply(strsplit(LastCommonAncestor_Name, ";"), tail, 1))

# Ahora usamos TaxonLevel para el heatmap
tabla_heat <- hmmer_annot %>%
  group_by(sample, TaxonLevel) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = TaxonLevel, values_from = count, values_fill = 0)

matriz <- as.matrix(tabla_heat[,-1])
rownames(matriz) <- tabla_heat$sample

matriz_lim <- pmin(matriz, 5)

pheatmap(matriz_lim,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         fontsize = 8,
         fontsize_row = 8,
         fontsize_col = 8,
         angle_col = 45,
         main = "Distribución taxonómica de VOGs contigs virales segun VIBRANT")

