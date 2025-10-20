library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(igraph)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(readr)
library(tibble)
library(stringr)

# -------------------------
# 1. Carpeta con archivos GTDB-Tk
# -------------------------
ruta_carpeta <- "/home/jesus/CH_22/galaxy_downloads/gtdb"
ruta_carpeta_egg <-"/home/jesus/CH_22/galaxy_downloads/eggnog Bins no virales, jesus1234"
ko_list <- read_delim("/home/jesus/CH_22/metagenoma/ko_list.txt", delim = "\t", col_names = TRUE)

archivos_tabular <- list.files(path = ruta_carpeta, pattern = "\\.tsv$", 
                               full.names = TRUE, recursive = TRUE)

archivos_tabular_egg <- list.files(path=ruta_carpeta_egg, pattern = "\\.tabular",
                                   full.names = TRUE, recursive = TRUE)


archivos_no_vacios <- archivos_tabular[file.size(archivos_tabular) > 0]

archivos_no_vacios_egg <- archivos_tabular_egg[file.size(archivos_tabular_egg) > 0]


# Crear nombres únicos para cada archivo usando carpeta + nombre de archivo
nombres_archivos <- sapply(archivos_no_vacios, function(x) {
  paste(basename(dirname(x)), tools::file_path_sans_ext(basename(x)), sep="_")
})

nombres_archivos_egg <- sapply(archivos_no_vacios_egg, function(x) {
  paste(basename(dirname(x)), tools::file_path_sans_ext(basename(x)), sep="_")
})

# Leer cada archivo como dataframe y guardarlo en una lista
lista_dataframes <- lapply(archivos_no_vacios, read.delim, header = TRUE, stringsAsFactors = FALSE)
names(lista_dataframes) <- nombres_archivos

lista_dataframes_egg <- lapply(archivos_no_vacios_egg, read.delim, header = TRUE, stringsAsFactors = FALSE)
names(lista_dataframes_egg) <- nombres_archivos_egg


# -------------------------
# 2. Revisar la lista
# -------------------------

length(lista_dataframes)  # número de archivos leídos
names(lista_dataframes)   # nobres de cada dataframe

length(lista_dataframes_egg)
names(lista_dataframes_egg)

lista_dataframes[[1]]   # primer archivo
lista_dataframes[["BM_F_summary"]]  # si ese fuera el nombre asignado

#funcion para extraer datos

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

# ------------------------------
# Extraer datos relevantes
# ------------------------------
columnas_por_df <- lapply(lista_dataframes_egg, names)

# Combinar todos los nombres en un solo vector
todas_las_columnas <- unique(unlist(columnas_por_df))

todas_las_columnas

ko_numero <- extraer_datos(lista_dataframes_egg,c("KEGG_ko","X.query"))

head(ko_numero)

# -------------------------
# 1. Filtrar dataframes bac120.summary
# -------------------------

df_bac120 <- lista_dataframes[grep("bac120\\.summary$", names(lista_dataframes))]

head(df_bac120)

df_ar53 <- lista_dataframes[grep("ar53\\.summary", names(lista_dataframes))]

head(df_ar53)
# -------------------------
# 2. Crear lista de tablas de abundancia por muestra
# -------------------------

df_total <- c(df_bac120, df_ar53)


taxo_list <- list()

for (muestra in names(df_total)) {
  df <- df_total[[muestra]]
  
  if ("classification" %in% colnames(df)) {
    # Extraer clase y orden
    clasificaciones <- sapply(strsplit(df$classification, ";"), function(x) {
      c_entry <- grep("^c__", x, value = TRUE)
      o_entry <- grep("^o__", x, value = TRUE)

      # Combinar en formato c__Clase_o__Orden
      if (length(c_entry) == 0) c_entry <- "c__Unclassified"
      if (length(o_entry) == 0) o_entry <- "o__Unclassified"
      paste0(c_entry, "_", o_entry)
    })
    
    # Eliminar sin clasificación
    clasificaciones <- clasificaciones[!grepl("Unclassified", clasificaciones)]
    
    if (length(clasificaciones) > 0) {
      counts <- as.data.frame(table(clasificaciones))
      colnames(counts) <- c("classification", muestra)
      taxo_list[[muestra]] <- counts
    }
  }
}

# -------------------------
# 3. Combinar todas las tablas en una matriz
# -------------------------
taxo_merged <- Reduce(function(x, y) merge(x, y, by="classification", all=TRUE), taxo_list)

# Reemplazar NAs por 0
taxo_merged[is.na(taxo_merged)] <- 0

# Limpieza de nombres de columna (muestras)
colnames(taxo_merged) <- gsub("^binning_", "", colnames(taxo_merged))          # elimina 'binning_'
colnames(taxo_merged) <- gsub("_gtdbtk", "", colnames(taxo_merged))            # elimina '_gtdbtk'
colnames(taxo_merged) <- gsub("\\.summary$", "", colnames(taxo_merged))  
rownames(taxo_merged) <- taxo_merged$classification
taxo_matrix <- as.matrix(taxo_merged[,-1])

# Convertir matriz a formato largo
taxo_long <- melt(taxo_matrix)
colnames(taxo_long) <- c("class", "sample", "count")
taxo_long

# Calcular porcentaje por muestra
taxo_long <- taxo_long %>%
  group_by(sample) %>%
  mutate(pct = count / sum(count) * 100)

taxo_long_filtrado<-taxo_long %>% filter(count > 0)

head(taxo_long_filtrado)
taxo_long_filtrado <- taxo_long_filtrado %>%
  filter(!grepl("^c___o__$", class)) %>%       # elimina c___o__ exacto
  filter(!grepl("_o__$", class))

#grafico de burbujas 

ggplot(taxo_long_filtrado, aes(x = sample, y = class, size = count, color = class)) +
  geom_point(alpha = 1) +
  scale_size_area(max_size = 15) +
  theme_classic() +  # elimina la cuadrícula de fondo
  labs(x = "Sample", y = "Class", size = "Abundance") +  # quitamos 'color = "Class"'
  guides(color = "none") +  # elimina la leyenda de color
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 25, face = "bold")
  ) +
  ggtitle("Taxonomic assignment to bins (GTDB-tk)")


# --- 1️⃣ Crear columna de bin base en taxonomía ---
taxo_long_filtrado <- taxo_long_filtrado %>%
  mutate(bin = str_extract(sample, "binning_[A-Z]{2}_[FS]"))  # Ej: binning_CQ_F

head(taxo_long_filtrado)

# --- 2️⃣ Crear columna de bin base en funciones ---
ko_numero <- ko_numero %>%
  mutate(bin = str_extract(muestra, "binning_[A-Z]{2}_[FS]"))  # idem

head(ko_numero)
# --- 3️⃣ Unir los dos data frames ---
ko_taxo <- left_join(ko_numero, taxo_long_filtrado, by = "bin")

# --- 4️⃣ Verificar ---
head(ko_taxo)


#------------------------
# Union por KO y grafica kegg pathways


####
#------------
# 1er grafico keggrest

library(KEGGREST)
library(pathview)
listDatabases()
library(dplyr)
library(tidyr)
library(RColorBrewer)

colnames(ko_numero)
head(ko_numero)

# --- 1. Contar KOs por muestra ---
ko_counts <- ko_numero %>%
  group_by(muestra, KEGG_ko) %>%
  summarise(Abundance = n(), .groups = "drop")

head(ko_counts)

# --- 2. Obtener todos los enlaces KO → pathway ---
all_links <- keggLink("pathway", "ko")
all_links_df <- data.frame(
  KO = sub("ko:", "", names(all_links)),
  Pathway = sub("path:", "", all_links),
  stringsAsFactors = FALSE
)

head(ko_counts$KEGG_ko) ##formato extraño no sirve para hacer la union con all_links

# --- Limpiar columna KEGG_ko ---

ko_counts_limpio <- ko_counts %>%
  # Separar múltiples KOs en filas distintas (por si hay "ko:K00887,ko:K00901")
  separate_rows(KEGG_ko, sep = ",") %>%
  # Quitar el prefijo "ko:"
  mutate(KEGG_ko = sub("^ko:", "", KEGG_ko)) %>%
  filter(KEGG_ko != "-")  # Eliminar los que no tienen KO válido

head(ko_counts_limpio) #se debe por ultimo renombrar la columna KEGG_ko  a KO 

ko_counts_limpio <- ko_counts_limpio %>%
  rename(KO = KEGG_ko)
 
head(ko_counts_limpio)  #listo 

# --- 3. Filtrar solo los KO presentes en tu dataset ---
ko_to_pathway <- subset(all_links_df, KO %in% ko_counts_limpio$KO)

head(ko_to_pathway)

# --- 4. Unir abundancias con pathways ---
ko_table <- merge(ko_counts_limpio, ko_to_pathway, by = "KO")

head(ko_table)  # ahora tienes KO, muestra, Abundance, Pathway

# --- 5. Sumar abundancia por pathway y muestra ---
pathway_summary <- ko_table %>%
  group_by(muestra, Pathway) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

head(pathway_summary)

# --- 6. Añadir nombres de pathways ---
path_names <- keggList("pathway")
path_names_df <- data.frame(
  Pathway = sub("path:", "", names(path_names)),
  Description = as.character(path_names),
  stringsAsFactors = FALSE
)

pathway_summary_named <- merge(pathway_summary, path_names_df, by = "Pathway")

head(pathway_summary_named)

pathways_filtered <- pathway_summary_named %>%
  filter(Description != "Metabolic pathways") 

pathways_filtered

###########################################################################

top_pathways <- pathways_filtered %>%
  group_by(muestra) %>%
  slice_max(order_by = Abundance, n = 1) %>%  # selecciona la de mayor abundancia
  ungroup()

top_pathways

# --- 1. Transformar los datos ---
matriz_df <- top_pathways %>%
  select(muestra, Abundance, Description) %>%
  pivot_wider(
    names_from = muestra,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  as.data.frame()   # data.frame clásico

# --- 2. Asignar Description como rownames ---
rownames(matriz_df) <- matriz_df$Description

# --- 3. Seleccionar solo columnas numéricas (las abundancias) ---
matriz <- as.matrix(matriz_df[, sapply(matriz_df, is.numeric)])

# --- 4. Heatmap ---
pheatmap(
  matriz,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  fontsize_row = 12,
  fontsize_col = 10,
  angle_col = 90,
  main = "Most commmon metabolic KEGG pathways per sample VIBRANT viruses",
  legend = TRUE
)
###Agregar al grafico de arriba la asignacion de hospedadores o el nombre del virus segun genomad o hmmer 

#---------------------
#Sin aplicar filtros
###

pathways_filtered

# --- 1. Transformar los datos ---

matriz_df <- pathways_filtered %>%
  select(muestra, Abundance, Description) %>%
  pivot_wider(
    names_from = muestra,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  as.data.frame()   # data.frame clásico

head(matriz_df)
# --- 2. Asignar Description como rownames ---
rownames(matriz_df) <- matriz_df$Description

# --- 3. Seleccionar solo columnas numéricas (las abundancias) ---
matriz <- as.matrix(matriz_df[, sapply(matriz_df, is.numeric)])

# --- 4. Heatmap ---
pheatmap(
  matriz,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = 90,
  main = "Most commmon metabolic KEGG pathways per sample VIBRANT viruses",
  legend = TRUE
)

#-------------------------
###Ahora buscaremos vias metabolicas mediante palabras clave 

# Lista de keywords de interés
keywords <- c("iron","carbon", "nitrate")

# Filtrar pathways con esas palabras clave (insensible a mayúsculas)
pathways_macro <- pathways_filtered %>%
  filter(grepl(paste(keywords, collapse = "|"), Description, ignore.case = TRUE))

# Revisar resultados
pathways_macro

# Agrupar y sumar abundancia por muestra y Description
macro_summary <- pathways_macro %>%
  group_by(muestra, Description) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = "drop")

macro_summary

matriz_df_2 <- macro_summary %>%
  select(muestra, Total_Abundance, Description) %>%
  pivot_wider(
    names_from = muestra,
    values_from = Total_Abundance,
    values_fill = 0
  ) %>%
  as.data.frame()   # data.frame clásico

# --- 2. Asignar Description como rownames ---
rownames(matriz_df_2) <- matriz_df_2$Description

# --- 3. Seleccionar solo columnas numéricas (las abundancias) ---
matriz <- as.matrix(matriz_df_2[, sapply(matriz_df_2, is.numeric)])

# --- 4. Heatmap ---
pheatmap(
  matriz,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = 90,
  main = "Most commmon metabolic KEGG pathways per sample VIBRANT viruses",
  legend = TRUE
)



#------------------------------------------------------------
####UNI los 2 data_frame en 1 segun los bins a los que correspondian 

colnames(ko_taxo)
head(ko_taxo)


# --- 1. Contar KOs por muestra ---
ko_counts_3 <- ko_taxo %>%
  group_by(X.query, KEGG_ko, bin) %>%
  summarise(Abundance = n(), class, .groups = "drop")

head(ko_counts_3)

# --- 2. Obtener todos los enlaces KO → pathway ---
all_links <- keggLink("pathway", "ko")
all_links_df <- data.frame(
  KO = sub("ko:", "", names(all_links)),
  Pathway = sub("path:", "", all_links),
  stringsAsFactors = FALSE
)

head(ko_counts_3$KEGG_ko) ##formato extraño no sirve para hacer la union con all_links

# --- Limpiar columna KEGG_ko ---

ko_counts_limpio_3 <- ko_counts_3 %>%
  # Separar múltiples KOs en filas distintas (por si hay "ko:K00887,ko:K00901")
  separate_rows(KEGG_ko, sep = ",") %>%
  # Quitar el prefijo "ko:"
  mutate(KEGG_ko = sub("^ko:", "", KEGG_ko)) %>%
  filter(KEGG_ko != "-")  # Eliminar los que no tienen KO válido

head(ko_counts_limpio_3) #se debe por ultimo renombrar la columna KEGG_ko  a KO 

ko_counts_limpio_3 <- ko_counts_limpio_3 %>%
  rename(KO = KEGG_ko)

head(ko_counts_limpio_3)  #listo 

# --- 3. Filtrar solo los KO presentes en tu dataset ---
ko_to_pathway_3 <- subset(all_links_df, KO %in% ko_counts_limpio_3$KO)

head(ko_to_pathway_3)

# --- 4. Unir abundancias con pathways ---
ko_table <- merge(ko_counts_limpio_3, ko_to_pathway_3, by = "KO")

head(ko_table)  # ahora tienes KO, muestra, Abundance, Pathway

# --- 5. Sumar abundancia por pathway y muestra ---
pathway_summary <- ko_table %>%
  group_by(bin, Pathway, X.query) %>%
  reframe(
    Abundance = sum(Abundance),
    class = unique(class)
  )

head(pathway_summary)

# --- 6. Añadir nombres de pathways ---
path_names <- keggList("pathway")
path_names_df <- data.frame(
  Pathway = sub("path:", "", names(path_names)),
  Description = as.character(path_names),
  stringsAsFactors = FALSE
)

pathway_summary_named <- merge(pathway_summary, path_names_df, by = "Pathway")

head(pathway_summary_named)

pathways_filtered <- pathway_summary_named %>%
  filter(Description != "Metabolic pathways") 

pathways_filtered

###########################################################################

top_pathways <- pathways_filtered %>%
  group_by(bin) %>%
  slice_max(order_by = Abundance, n = 1) %>%  # selecciona la de mayor abundancia
  ungroup()

head(top_pathways)

###Hacer heatmap con anotacion taxonomica de los bins, mostrar bins como grupos en colores, nivel clases o sp de ser posible por muestra 

heat_data <- top_pathways %>%
  group_by(Description, bin, class,X.query) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

head(heat_data)
library(reshape2)

heat_mat_df <- dcast(
  heat_data,
  Description ~ bin,
  value.var = "Abundance",
  fun.aggregate = sum,
  fill = 0
)

# Fijar nombres de fila y convertir a matriz
rownames(heat_mat_df) <- heat_mat_df$Description
heat_mat <- as.matrix(heat_mat_df[, -1, drop = FALSE])

# --- Heatmap básico ---
pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Vías metabólicas por bin"
)

# crear anotación (por bin)
ann_col <- heat_data %>%
  group_by(bin) %>%
  summarise(class = first(na.omit(class)), .groups = "drop") %>%
  as.data.frame()

rownames(ann_col) <- ann_col$bin
ann_col <- ann_col[, "class", drop = FALSE]

# (opcional) definir colores para cada clase
clases_unicas <- unique(na.omit(ann_col$class))
paleta <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(clases_unicas)),
  clases_unicas
)
ann_colors <- list(class = paleta)

pheatmap(
  heat_mat,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "row",
  fontsize_row = 12,
  fontsize_col = 12,
  main = "Vías metabólicas por bin (agrupadas por clase)"
)





