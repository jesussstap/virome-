library(dplyr)
library(tidyr)
library(pheatmap)
library(readr)
library(tools)
# --- 1. Leer archivos HMMER ---

ruta <- "/home/jesus/CH_22/galaxy_downloads/iphop/todos_virus/todos_los_virus"
archivos <- list.files(ruta, pattern = "\\.csv$", full.names = TRUE)
archivos


leer_tabla <- function(f) {
  df <- suppressMessages(read_delim(f, delim = NULL, comment = "#"))  # autodetecta
  if (nrow(df) == 0) return(NULL)
  df$sample <- file_path_sans_ext(basename(f))
  return(df)
}

  
lista_res <- lapply(archivos, leer_tabla)

extraer_datos <- function(lista_df, columnas) {
  # validación simple
  if (missing(lista_df) || length(lista_df) == 0) return(NULL)
  if (missing(columnas)) stop("'columnas' debe indicarse (vector de nombres de columna).")
  columnas <- as.character(columnas)
  
  # asegurar nombres para las muestras
  if (is.null(names(lista_df))) names(lista_df) <- paste0("sample_", seq_along(lista_df))
  
  resultados <- list()
  
  for (i in seq_along(lista_df)) {
    df <- lista_df[[i]]
    nombre_muestra <- names(lista_df)[i]
    
    # columnas de interés que están en este dataframe
    cols_presentes <- intersect(columnas, names(df))
    if (length(cols_presentes) == 0) next    # nada que extraer de este df
    
    tmp <- df[, cols_presentes, drop = FALSE]
    
    # si no hay filas (dataframe vacío), saltar
    if (nrow(tmp) == 0) next
    
    # Añadir columnas solicitadas que faltan (con NA), para homogeneizar
    faltantes <- setdiff(columnas, names(tmp))
    if (length(faltantes) > 0) {
      for (mc in faltantes) tmp[[mc]] <- NA
    }
    
    # Reordenar columnas según el orden pedido
    tmp <- tmp[, columnas, drop = FALSE]
    
    # Añadir columna 'muestra' con el nombre repetido por fila
    tmp$muestra <- rep(nombre_muestra, nrow(tmp))
    
    resultados[[length(resultados) + 1]] <- tmp
  }
  
  if (length(resultados) == 0) return(NULL)
  
  # unir todo en un solo data.frame
  dplyr::bind_rows(resultados)
}


hospedador<-extraer_datos(lista_res, c("Virus", "Host taxonomy","Confidence score"))

#filtrar hospedadores con mayor confidence score, 1 por node, existen nodes con mas de 1 hospedador

confianza_mas_alta_x_node <- hospedador %>%
  group_by(Virus) %>%
  slice_max(order_by = `Confidence score`, n = 1, with_ties = FALSE) %>%  # elegir la fila de mayor confianza
  ungroup() %>%
  select(muestra, Virus, `Host taxonomy`, `Confidence score`)

head(confianza_mas_alta_x_node)

#para ver si efectivamente nos quedamos con solo 1 host por virus, deberia dar 0

print (confianza_mas_alta_x_node %>%
  count(Virus) %>%
  filter(n > 1))

hospedador_2 <- confianza_mas_alta_x_node %>%
  filter(!is.na(`Host taxonomy`)) %>%
  mutate(TaxonLevel = sapply(strsplit(`Host taxonomy`, ";"), tail, 1))

head(hospedador_2$TaxonLevel)


#Datos de virus 

##--------------------------###
# Hacer hmmer de nuevo, agregar numero de node 
##--------------------------###

ruta_hmmer <- "/home/jesus/CH_22/metagenoma/resultados_hmmer"
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

lista_res_vir <- lapply(archivos, leer_hmmer)
hmmer_res <- do.call(rbind, lista_res_vir[!sapply(lista_res, is.null)])


# --- 2. Leer tabla de anotaciones ---
vog_annot <- read.delim("/home/jesus/CH_22/metagenoma/vog_hmm/vog.lca.tsv",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

colnames(vog_annot)[1] <- "target_name"

# --- 3. Unir resultados de HMMER con anotaciones ---
hmmer_annot <- hmmer_res %>%
  left_join(vog_annot, by = "target_name")

# --- 4. Preparar tabla para heatmap usando taxonomía ---
tax_col <- "LastCommonAncestor_Name"  # categoría para el heatmap

# Extraer solo el nivel más bajo de LastCommonAncestor_Name
hmmer_annot <- hmmer_annot %>%
  filter(!is.na(LastCommonAncestor_Name)) %>%
  mutate(TaxonLevel_viral = sapply(strsplit(LastCommonAncestor_Name, ";"), tail, 1))


# Ahora usamos TaxonLevel para el heatmap
tabla_heat <- hospedador_2 %>%
  group_by(muestra, TaxonLevel) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = TaxonLevel, values_from = count, values_fill = 0)

matriz <- as.matrix(tabla_heat[,-1])
rownames(matriz) <- tabla_heat$muestra

matriz_lim <- pmin(matriz, 5)

pheatmap(matriz_lim,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         fontsize = 8,
         fontsize_row = 8,
         fontsize_col = 8,
         angle_col = 45,
         main = "Heatmap hospedadores virus")




library(dplyr)
library(tidyr)
library(igraph)
library(proxy)

# -------------------
# 1. Matriz de abundancias hospedador vs muestra
# -------------------

hosp_matrix <- hospedador_2 %>%
  group_by(Virus, TaxonLevel) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Virus, values_from = count, values_fill = 0)

# Filas = hospedadores (TaxonLevel), columnas = muestras

matriz_binaria <- (hosp_matrix[, -1] > 0) * 1
rownames(matriz_binaria) <- hosp_matrix$TaxonLevel

# -------------------
# 2. Similaridad entre muestras (transponemos para comparar muestras)
# -------------------

similaridad <- proxy::simil(t(matriz_binaria), method = "bray")
adj <- as.matrix(similaridad)

# Opcional: aplicar un umbral (ejemplo, >0.3)
umbral <- 0.2
adj[adj < umbral] <- 0

# -------------------
# 3. Construcción del grafo
# -------------------
g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
V(g)$name <- colnames(matriz_binaria)
V(g)$label <- V(g)$name

# Colores por grupo de muestra (ej: prefijo antes del "_")
grupos <- sub("_.*", "", V(g)$name)
paleta <- rainbow(length(unique(grupos)))
V(g)$color <- paleta[as.factor(grupos)]

# Layout
lay <- layout_with_fr(g)

# -------------------
# 4. Plot
# -------------------
plot(
  g,
  layout = lay,
  vertex.size = 20,
  vertex.label.cex = 0.8,
  vertex.label.color = "black",
  edge.width = E(g)$weight * 2,
  main = "Red de similaridad entre muestras (hospedadores)"
)

legend("topleft",
       legend = unique(grupos),
       col = paleta,
       pch = 19,
       pt.cex = 1.5,
       bty = "n",
       title = "Grupos de muestras")

# Si quieres explorarlo en ventana interactiva
tkplot(g, canvas.width = 500, canvas.height = 500)








#Unir hmmeranot con hospedador_2 por numero de contig

##cambie el nombre de los contigs, prokka no leia los contigs, necestio que los  
## node esten asociados a la secuencia de aa resultante.


#Otra manera puede ser simplemente extraer el nombre de los virus con hospedadores reconocidos y luego hacer hmmer a parte 

library(Biostrings)
library(stringr)

# Ruta donde están los fasta en subcarpetas
ruta_fasta <- "/home/jesus/CH_22/galaxy_downloads/vibrant jesus.tapia"

# Buscar todos los archivos .fasta en subcarpetas
archivos_fasta <- list.files(ruta_fasta, pattern = "\\.fasta$", 
                             full.names = TRUE, recursive = TRUE)

length(archivos_fasta)  # cuántos encontró

# Lista con todos los fasta leídos
lista_fasta <- lapply(archivos_fasta, readDNAStringSet)

# Poner nombre a cada objeto según el archivo
names(lista_fasta) <- basename(archivos_fasta)

# Vector con los nombres de virus que quieres extraer
nombre_virus <- unique(hospedador_2$Virus)

# Buscar secuencias cuyo nombre contenga alguno de los virus seleccionados
extraer_fasta <- function(dna, virus_names) {
  idx <- grepl(paste(virus_names, collapse = "|"), names(dna))
  dna[idx]
}

# Aplicar a cada archivo fasta
fasta_filtrados <- lapply(lista_fasta, extraer_fasta, virus_names = nombre_virus)

fasta_no_vacios <- fasta_filtrados[sapply(fasta_filtrados, length) > 0]

# Concatenar todos los DNAStringSet en uno solo
fasta_final <- unlist(DNAStringSetList(fasta_no_vacioss))

# Confirmar clase
print(class(fasta_final))     # ahora debería ser "DNAStringSet"
print(length(fasta_final))    # número de secuencias totales

# Guardar a archivo FASTA
output_file <- "/home/jesus/CH_22/virus_con_hospedador.fasta"
writeXStringSet(fasta_final, filepath = output_file)



#### Lo que hay que hacer es buscar los nodes (contig viral) y ver a que muestra corresponde para saber donde encontramos esos hospedadores

#hacer red de similitud entre muestras segun los virus encontrados y su hospedador 
#aplicar distancia de jaccard 





# Librerías necesarias
library(igraph)
library(Hmisc)   # Para calcular correlaciones y p-valores
library(reshape2)

# 1. Cargar tu tabla de abundancias (filas = virus, columnas = muestras)
# Ejemplo: abundancias virales
# Cada fila: virus/contig, cada columna: muestra
abund <- read.csv("virus_abundancias.csv", row.names = 1)

# 2. Calcular correlaciones de co-ocurrencia
res <- rcorr(as.matrix(t(abund)), type = "spearman")  # o "pearson"

# Extraer correlaciones y p-valores
corr <- res$r
pval <- res$P

# 3. Filtrar correlaciones significativas
umbral_cor <- 0.6   # solo correlaciones fuertes (>0.6)
umbral_p <- 0.05    # significancia estadística

# Matriz filtrada
corr_filtrada <- corr
corr_filtrada[abs(corr) < umbral_cor | pval > umbral_p] <- 0

# 4. Crear red
g <- graph_from_adjacency_matrix(corr_filtrada, mode = "undirected", weighted = TRUE, diag = FALSE)

# 5. Visualizar red
plot(g,
     vertex.size = 8,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.color = "skyblue",
     edge.width = E(g)$weight * 2,
     edge.color = ifelse(E(g)$weight > 0, "darkgreen", "red")
)
