# ============================
#   Vibrant - Análisis gráfico
# ============================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")  # instala BiocManager si no lo tienes
}
BiocManager::install("ggraph")
BiocManager::install("Biostrings")
install.packages("proxy","igraph","tridytext","ggraph")
library(Biostrings)
library(proxy)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(treemapify)
library(tidytext)
library(ggpubr)
library(igraph)
library(vegan)
library(ggraph)


# Ruta con archivos de salida VIBRANT
ruta_carpeta <- "/home/jesus/CH_22/galaxy_downloads/vibrant_jesus.tapia"

# Obtener archivos .tabular no vacíos
archivos_tabular <- list.files(path = ruta_carpeta, pattern = "\\.tabular$", 
                               full.names = TRUE, recursive = TRUE)

archivos_no_vacios <- archivos_tabular[file.size(archivos_tabular) > 0]

if (length(archivos_no_vacios) == 0) {
  stop("No se encontraron archivos con datos en la carpeta especificada.")
}

# Extraer nombre de la carpeta (muestra)
nombres_muestra <- basename(dirname(archivos_no_vacios))

# Leer archivos como lista de dataframes
lista_dataframes <- lapply(archivos_no_vacios, read.delim, header = TRUE)
names(lista_dataframes) <- nombres_muestra   # ahora los nombres son las carpetas

# ------------------------------
# Función para extraer columna con ID de muestra
# ------------------------------
#Para extraer 1 solo dato 
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
# Extraer nombres de columnas de cada dataframe
columnas_por_df <- lapply(lista_dataframes, names)

# Combinar todos los nombres en un solo vector
todas_las_columnas <- unique(unlist(columnas_por_df))

# Mostrar el resultado
print(todas_las_columnas)

sca<-extraer_datos(lista_dataframes, "scaffold")
contig_calidad<-extraer_datos(lista_dataframes,c("scaffold","Quality"))
metadata_genoma_completo<- extraer_datos(lista_dataframes,c("protein",                        "scaffold",                          "KO"     ,                            
                                                             "AMG",                                "KO.name",                             "KO.evalue"  ,                        
                                                             "KO.score",                            "KO.v.score",                          "Pfam"       ,                        
                                                             "Pfam.name",                           "Pfam.evalue",                         "Pfam.score"  ,                       
                                                             "Pfam.v.score",                        "VOG",                                 "VOG.name"     ,                      
                                                             "VOG.evalue"   ,                       "VOG.score",                           "VOG.v.score"   ,                     
                                                             "fragment"      ,                      "protein.start",                       "protein.stop"   ,                    
                                                             "protein.length" ,                     "nucleotide.start",                    "nucleotide.stop" ,                   
                                                             "nucleotide.length",                   "total.genes",                         "all.KEGG"         ,                  
                                                             "KEGG.v.score"      ,                  "all.Pfam",                            "all.VOG"           ,                 
                                                             "KEGG.int.rep"      ,                 "KEGG.zero",                           "Pfam.int.rep"        ,         
                                                             "Pfam.zero"          ,                 "VOG.redoxin",                         "VOG.rec.tran"         ,              
                                                             "VOG.int"             ,                "VOG.RnR",                             "VOG.DNA"               ,             
                                                             "KEGG.restriction.check",              "KEGG.toxin.check",                    "VOG.special"            ,            
                                                             "annotation.check"       ,             "p_v.check",                           "p_k.check"               ,           
                                                             "k_v.check"               ,            "k.check",                             "p.check"                  ,          
                                                             "v.check"                  ,           "h.check" ,                            "type"                      ,         
                                                             "Quality"                   ,          "AMG.count",                           "AMG.KO"   ,                          
                                                             "AMG.KO.name"                ,         "prediction",                           
                                                             "lytic"                       ,        "complete.circular"))
metadata_genoma_completo
prediccion_df <- extraer_datos(lista_dataframes, "prediction")
tipo_virus_df <- extraer_datos(lista_dataframes, "type")
vog_name_df   <- extraer_datos(lista_dataframes, "VOG.name")
calidad_df    <- extraer_datos(lista_dataframes, "Quality")
AMG_df   <- extraer_datos(lista_dataframes, "AMG")
genoma_circular <- extraer_datos(lista_dataframes, "Quality")
Ko_name<- extraer_datos(lista_dataframes, "KO.name" )
Pfam_scan<- extraer_datos(lista_dataframes, "Pfam.name")
contig<- extraer_datos(lista_dataframes, "protein")#protein es el nombre del  
amgs <- extraer_datos(lista_dataframes,c("AMG.KO.name"))

write.table(amgs,
            file = "/home/jesus/CH_22/metagenoma/AMG_VIBRANT/amgs.tabular",  # ruta y nombre del archivo
            sep = "\t",       # separador tabulación
            row.names = FALSE, # no guardar nombres de fila
            quote = FALSE) 


# contig, aqui fallo el programa anterior al asignar el nombre a la tabla

circular_ids <- metadata_genoma_completo %>%
  filter(Quality == "complete circular") %>%
  select(muestra,"protein",                        "scaffold",                          "KO"     ,                            
         "AMG",                                "KO.name",                             "KO.evalue"  ,                        
         "KO.score",                            "KO.v.score",                          "Pfam"       ,                        
         "Pfam.name",                           "Pfam.evalue",                         "Pfam.score"  ,                       
         "Pfam.v.score",                        "VOG",                                 "VOG.name"     ,                      
         "VOG.evalue"   ,                       "VOG.score",                           "VOG.v.score"   ,                     
         "fragment"      ,                      "protein.start",                       "protein.stop"   ,                    
         "protein.length" ,                     "nucleotide.start",                    "nucleotide.stop" ,                   
         "nucleotide.length",                   "total.genes",                         "all.KEGG"         ,                  
         "KEGG.v.score"      ,                  "all.Pfam",                            "all.VOG"           ,                 
         "KEGG.int.rep"      ,                 "KEGG.zero",                           "Pfam.int.rep"        ,         
         "Pfam.zero"          ,                 "VOG.redoxin",                         "VOG.rec.tran"         ,              
         "VOG.int"             ,                "VOG.RnR",                             "VOG.DNA"               ,             
         "KEGG.restriction.check",              "KEGG.toxin.check",                    "VOG.special"            ,            
         "annotation.check"       ,             "p_v.check",                           "p_k.check"               ,           
         "k_v.check"               ,            "k.check",                             "p.check"                  ,          
         "v.check"                  ,           "h.check" ,                            "type"                      ,         
         "Quality"                   ,          "AMG.count",                           "AMG.KO"   ,                          
         "AMG.KO.name"                ,         "prediction",                           
         "lytic"                       ,        "complete.circular")   # reemplaza ID por la columna que tenga el identificador de tu genoma/bin

circular_ids
ko_numero<- extraer_datos(lista_dataframes,"KO")
ko_numero


high_q<- metadata_genoma_completo %>%
  filter(Quality == "high quality draft") %>%
  select(muestra, "protein",                        "scaffold",                          "KO"     ,                            
"AMG",                                "KO.name",                             "KO.evalue"  ,                        
"KO.score",                            "KO.v.score",                          "Pfam"       ,                        
"Pfam.name",                           "Pfam.evalue",                         "Pfam.score"  ,                       
"Pfam.v.score",                        "VOG",                                 "VOG.name"     ,                      
"VOG.evalue"   ,                       "VOG.score",                           "VOG.v.score"   ,                     
"fragment"      ,                      "protein.start",                       "protein.stop"   ,                    
"protein.length" ,                     "nucleotide.start",                    "nucleotide.stop" ,                   
"nucleotide.length",                   "total.genes",                         "all.KEGG"         ,                  
"KEGG.v.score"      ,                  "all.Pfam",                            "all.VOG"           ,                 
"KEGG.int.rep"      ,                 "KEGG.zero",                           "Pfam.int.rep"        ,         
"Pfam.zero"          ,                 "VOG.redoxin",                         "VOG.rec.tran"         ,              
"VOG.int"             ,                "VOG.RnR",                             "VOG.DNA"               ,             
"KEGG.restriction.check",              "KEGG.toxin.check",                    "VOG.special"            ,            
"annotation.check"       ,             "p_v.check",                           "p_k.check"               ,           
"k_v.check"               ,            "k.check",                             "p.check"                  ,          
"v.check"                  ,           "h.check" ,                            "type"                      ,         
"Quality"                   ,          "AMG.count",                           "AMG.KO"   ,                          
"AMG.KO.name"                ,         "prediction",                           
"lytic"                       ,        "complete.circular")

write.table(high_q,
            file = "/home/jesus/CH_22/alta_calidad/alta_calidad.tabular",  # ruta y nombre del archivo
            sep = "\t",       # separador tabulación
            row.names = FALSE, # no guardar nombres de fila
            quote = FALSE) 

media_q<-metadata_genoma_completo %>%
    filter(Quality == "medium quality draft") %>%
    select(muestra, "protein",                        "scaffold",                          "KO"     ,                            
           "AMG",                                "KO.name",                             "KO.evalue"  ,                        
           "KO.score",                            "KO.v.score",                          "Pfam"       ,                        
           "Pfam.name",                           "Pfam.evalue",                         "Pfam.score"  ,                       
           "Pfam.v.score",                        "VOG",                                 "VOG.name"     ,                      
           "VOG.evalue"   ,                       "VOG.score",                           "VOG.v.score"   ,                     
           "fragment"      ,                      "protein.start",                       "protein.stop"   ,                    
           "protein.length" ,                     "nucleotide.start",                    "nucleotide.stop" ,                   
           "nucleotide.length",                   "total.genes",                         "all.KEGG"         ,                  
           "KEGG.v.score"      ,                  "all.Pfam",                            "all.VOG"           ,                 
           "KEGG.int.rep"      ,                 "KEGG.zero",                           "Pfam.int.rep"        ,         
           "Pfam.zero"          ,                 "VOG.redoxin",                         "VOG.rec.tran"         ,              
           "VOG.int"             ,                "VOG.RnR",                             "VOG.DNA"               ,             
           "KEGG.restriction.check",              "KEGG.toxin.check",                    "VOG.special"            ,            
           "annotation.check"       ,             "p_v.check",                           "p_k.check"               ,           
           "k_v.check"               ,            "k.check",                             "p.check"                  ,          
           "v.check"                  ,           "h.check" ,                            "type"                      ,         
           "Quality"                   ,          "AMG.count",                           "AMG.KO"   ,                          
           "AMG.KO.name"                ,         "prediction",                           
           "lytic"                       ,        "complete.circular")


write.table(media_q,
            file = "/home/jesus/CH_22/media_calidad/media_calidad.tabular",  # ruta y nombre del archivo
            sep = "\t",       # separador tabulación
            row.names = FALSE, # no guardar nombres de fila
            quote = FALSE)

#sacar solo virus lisogenicos
lisogenicos<-metadata_genoma_completo %>%
  filter( type == "lysogenic") %>%
  select(muestra, "protein",                        "scaffold",                          "KO"     ,                            
         "AMG",                                "KO.name",                             "KO.evalue"  ,                        
         "KO.score",                            "KO.v.score",                          "Pfam"       ,                        
         "Pfam.name",                           "Pfam.evalue",                         "Pfam.score"  ,                       
         "Pfam.v.score",                        "VOG",                                 "VOG.name"     ,                      
         "VOG.evalue"   ,                       "VOG.score",                           "VOG.v.score"   ,                     
         "fragment"      ,                      "protein.start",                       "protein.stop"   ,                    
         "protein.length" ,                     "nucleotide.start",                    "nucleotide.stop" ,                   
         "nucleotide.length",                   "total.genes",                         "all.KEGG"         ,                  
         "KEGG.v.score"      ,                  "all.Pfam",                            "all.VOG"           ,                 
         "KEGG.int.rep"      ,                 "KEGG.zero",                           "Pfam.int.rep"        ,         
         "Pfam.zero"          ,                 "VOG.redoxin",                         "VOG.rec.tran"         ,              
         "VOG.int"             ,                "VOG.RnR",                             "VOG.DNA"               ,             
         "KEGG.restriction.check",              "KEGG.toxin.check",                    "VOG.special"            ,            
         "annotation.check"       ,             "p_v.check",                           "p_k.check"               ,           
         "k_v.check"               ,            "k.check",                             "p.check"                  ,          
         "v.check"                  ,           "h.check" , 
         "Quality"                   ,          "AMG.count",                           "AMG.KO"   ,                          
         "AMG.KO.name"                ,         "prediction",                           
         "lytic"                       ,        "complete.circular")


write.table(lisogenicos,
            file = "/home/jesus/CH_22/lisogenicos/lisogenicos.tabular",  # ruta y nombre del archivo
            sep = "\t",       # separador tabulación
            row.names = FALSE, # no guardar nombres de fila
            quote = FALSE)

####Para extraer plasmidos
plasmidos<- metadata_genoma_completo %>%
  filter(prediction == "plasmid") %>%
  select(muestra,"scaffold","prediction")


plasmidos
#Guardar carpeta con nombre de plasmidos y su muestra asociada 


write.table(plasmidos,
            file = "/home/jesus/CH_22/plasmidos_CH22/plasmidos.tabular",  # ruta y nombre del archivo
            sep = "\t",       # separador tabulación
            row.names = FALSE, # no guardar nombres de fila
            quote = FALSE) 

#Sacar genomas de alta calidad 



plasmidos
#Guardar carpeta con nombre de plasmidos y su muestra asociada 


write.table(plasmidos,
            file = "/home/jesus/CH_22/plasmidos_CH22/plasmidos.tabular",  # ruta y nombre del archivo
            sep = "\t",       # separador tabulación
            row.names = FALSE, # no guardar nombres de fila
            quote = FALSE) 

#con esto tenemos los nombres de las secuencias reconocidas como plasmidos, 
#ahora debemos obtener la secuencia.fasta correspondiente a cada plasmido,
#para esto debemos tener el archivo fasta original de entrada a vibrant, osea el 
#ensamble


plasmidos_CH22 <- data.frame(scaffold = c(plasmidos$scaffold))
plasmidos_CH22
muestra<- data.frame(muestra=c(plasmidos$muestra))

length(plasmidos_CH22$scaffold)
# --- 2. Carpeta con los FASTA ---
ruta_fasta <- "/home/jesus/metagenoma/ensamble_prueba/"

archivos_fasta <- list.files(
  path = ruta_fasta, 
  pattern = "\\.fasta$", 
  full.names = TRUE, 
  recursive = TRUE
)
archivos_fasta
######  Guardar plasmidos como archivo tabular en carpeta X ########## 

nombres_muestra <- tools::file_path_sans_ext(basename(archivos_fasta))
names(archivos_fasta) <- nombres_muestra

# --- 3. Función para extraer y guardar contigs individualmente ---
extraer_y_guardar_contigs <- function(muestra, plasmidos_CH22, archivos_fasta, ruta_salida) {
  if (!muestra %in% names(archivos_fasta)) {
    warning(paste("No se encontró archivo FASTA para muestra:", muestra))
    return(NULL)
  }
  
  fasta_file <- archivos_fasta[[muestra]]
  seqs <- readDNAStringSet(fasta_file)
  
  # Filtrar y guardar cada contig en un archivo separado
  for (scaf in scaffolds) {
    if (scaf %in% names(seqs)) {
      seq_filtrada <- seqs[scaf]
      nombre_archivo <- paste0(ruta_salida, "/", muestra, "_", scaf, ".fasta")
      writeXStringSet(seq_filtrada, filepath = nombre_archivo)
    }
  }
}

# --- 4. Crear carpeta de salida ---
ruta_salida <- "/home/jesus/CH_22/plasmidos_CH22"
if (!dir.exists(ruta_salida)) dir.create(ruta_salida)

# --- 5. Aplicar función a todas las muestras ---
for (m in unique(plasmidos$muestra)) {
  scaffolds <- plasmidos %>%
    filter(muestra == m) %>%
    pull(scaffold)
  
  extraer_y_guardar_contigs(m, scaffolds, archivos_fasta, ruta_salida)
}

#para ver si se crearon los archivos fasta de los plasmidos
length(plasmidos$scaffold)
list.files(
  path = ruta_salida, 
  pattern = "\\.fasta$", 
  full.names = TRUE, 
  recursive = TRUE
)
#Lo anterior indica que existen solo 76 secuencias desde las cuales se encontraron 
#plasmidos, pero como se observa al lista plasmidos_CH22$scaffold son 306 los fragmentos
#estas son secuencias y derivados de esas secuencias, por eso se infla el numero 
#trabajare con las 76 


####-------------------REVISAR ESTO
#   AMG_name_df <- extraer_datos(lista_dataframes, "AMG KO name") ##Se debe ingresar la tabla donde se encuentra este amg ko name 
####------------------
# ------------------------------
# GRÁFICO 1: Predicción (Virus vs Organismo)


prediccion_df <- extraer_datos(lista_dataframes, "prediction")

if (!is.null(prediccion_df)) {
  prediccion_clean <- prediccion_df %>%
    filter(!is.na(prediction)) %>%
    count(muestra, prediction, sort = TRUE)

  
  ggplot(prediccion_clean, aes(x = prediction, y = n, fill = muestra)) +
    geom_col(position = "dodge") +
    labs(title = "Viruses vs Organisms (VIBRANT)",
         x = "Prediction", y = "Contigs") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
} 

# ------------------------------
# GRÁFICO 2: Tipo de virus (lítico vs lisogénico)
# ------------------------------
if (!is.null(tipo_virus_df)) {
  tipo_clean <- tipo_virus_df %>%
    filter(!is.na(type)) %>%
    count(muestra, type, sort = TRUE) %>%
    filter(type %in% c("lytic", "lysogenic", "temperate"))
  
  if (nrow(tipo_clean) > 0) {
    ggplot(tipo_clean, aes(x = type, y = n, fill = muestra)) +
      geom_col(position = "dodge") +
      labs(title = "Viral type (VIBRANT)",
           x = "Type of viral cycle", y = "Contigs") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }
}

library("ggplot2")
library("ggrepel")
library("patchwork")
library(ggfortify)
library(dplyr)
library(tidyr)
library(readxl)


Datos_geoquimica_limpios_para_R <- read_excel("Datos geoquimica limpios para R.xlsx")
View(Datos_geoquimica_limpios_para_R)   
Datos<-Datos_geoquimica_limpios_para_R
Lugar<-Datos$muestra

tipo_clean

data_virus_geo <- bind_rows(Datos,tipo_clean)

colnames(data_virus_geo)

library(ggplot2)



ggplot(data_virus_geo, aes(x = pH , y = n, color = type)) +
  geom_point(size = 3, alpha = 0.7) +  # puntos más grandes y con transparencia
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") + # línea de regresión
  labs(
    title = "Relation betwen viral life cycle and pH (VIBRANT)",
    x = "pH",
    y = "Number of contigs",
    color = "Viral type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

library(lme4)

modelo <- lmer(abundancia ~ temp + pH + tipo_viral + (1|muestra), data = datos_largos)

summary(modelo)

# ------------------------------
# GRÁFICO 3: Top 20 proteínas VOG más abundantes
# ------------------------------

library(dplyr)
library(tidyr)
library(proxy)
library(igraph)
library(ggplot2)

# -------------------
# Limpiar y seleccionar top VOGs
# -------------------
vog_clean <- vog_name_df %>%
  mutate(VOG.short = ifelse(nchar(VOG.name) > 40,
                            paste0(substr(VOG.name, 1, 37), "..."),
                            VOG.name))

if (!is.null(vog_clean)) {
  vog_clean <- vog_clean %>%
    filter(!is.na(VOG.name), VOG.name != "",
           !grepl("hypothetical", VOG.name, ignore.case = TRUE)) %>%
    count(muestra, VOG.name, sort = TRUE) %>%
    group_by(muestra) %>%
    slice_max(n, n = 10)
  
  # Grafico de barras top VOGs
  if (nrow(vog_clean) > 0) {
    ggplot(vog_clean, aes(x = reorder(VOG.name, n), y = n, fill = muestra)) +
      geom_col(position = "dodge") +
      coord_flip() +
      labs(title = "Top 10 proteínas VOG más abundantes por muestra",
           x = "Proteína VOG", y = "Frecuencia") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.y = element_text(size = 7))
  }
}

# -------------------
# Crear matriz binaria VOG vs Muestra
# -------------------
VOG_matrix <- vog_clean %>%
  group_by(VOG.name) %>%
  summarise(total = sum(n)) %>%
  filter(total >= 1) %>%
  inner_join(vog_clean, by = "VOG.name") %>%
  pivot_wider(names_from = muestra, values_from = n, values_fill = 0)

# columnas 3 en adelante son muestras
matriz_binaria <- (VOG_matrix[, 3:ncol(VOG_matrix)] > 0) * 1
rownames(matriz_binaria) <- VOG_matrix$VOG.name

# -------------------
# Similaridad entre muestras (transponer para que columnas sean nodos)
# -------------------
similaridad <- proxy::simil(t(matriz_binaria), method = "bray")
adj <- as.matrix(similaridad)
umbral <- 0
adj[adj < umbral] <- 0

# -------------------
# Grafo: nodos = muestras
# -------------------
g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
V(g)$name <- colnames(matriz_binaria)  # nombres de nodos = muestras
V(g)$label <- V(g)$name

# Colores por grupo (prefijo de la muestra)
grupos <- sub("_.*", "", V(g)$name)
paleta <- rainbow(length(unique(grupos)))
V(g)$color <- paleta[as.factor(grupos)]

# Layout más informativo
lay <- layout_with_fr(g)

# -------------------
# Plot
# -------------------
plot(
  g,
  layout = lay,
  vertex.size = 20,
  vertex.label.cex = 0.8,
  vertex.label.color = "black",
  edge.width = E(g)$weight * 2,
  main = "Red de similaridad entre muestras (VOG checkV)"
)

legend("topleft",
       legend = unique(grupos),
       col = paleta,
       pch = 19,
       pt.cex = 1.5,
       bty = "n",
       title = "Grupos de muestras")

tkplot(g, canvas.width = 450, canvas.height = 450)


# ------------------------------
# GRÁFICO 5: Calidad de VIBRANT
# ------------------------------
if (!is.null(calidad_df)) {
  calidad_clean <- calidad_df %>%
    filter(!is.na(Quality)) %>%
    count(muestra, Quality, sort = TRUE)
  
  ggplot(calidad_clean, aes(x = Quality, y = n, fill = muestra)) +
    geom_col(position = "dodge") +
    labs(title = "Calidad de predicciones VIBRANT",
         x = "Calidad", y = "Número de scaffolds") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

#genomas circulares, se asume completos con todos sus metadatos
circular <- calidad_clean %>%
  filter(Quality == "complete circular") %>%
  count(muestra, Quality, sort = TRUE)

circular_ids <- contig_calidad %>%
  filter(Quality == "complete circular") %>%
  select(muestra, "scaffold")   # reemplaza ID por la columna que tenga el identificador de tu genoma/bin

circular_ids

#Grafico 6 AMG 

if (!is.null(AMG_df)) {
  AMG_clean <- AMG_df %>%
    filter(!is.na(AMG), AMG != "", AMG != "none") %>%  # <-- filtramos contigs sin AMGs
    count(muestra, AMG, sort = TRUE)
  
  ggplot(AMG_clean, aes(x = AMG, y = n, fill = muestra)) +
    geom_col(position = "dodge") +
    labs(title = "Distribución de AMGs por muestra",
         x = "Muestra", y = "Número de scaffolds") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

#Grafico 7 Heatmap Pfamscan Vibrant 

Pfam_scan

if (!is.null(Pfam_scan)) {
  Pfam_clean <- Pfam_scan %>%
    filter(!is.na(`Pfam.name`), `Pfam.name` != "") %>%
    count(muestra, `Pfam.name`, sort = TRUE) %>%
    group_by(muestra) %>%
    slice_max(n, n = 10)
  
  Pfam_matrix <- Pfam_clean %>%
    group_by(Pfam.name) %>%
    summarise(total = sum(n)) %>%
    filter(total >= 10) %>%      
    inner_join(Pfam_clean, by = "Pfam.name") %>%
    pivot_wider(names_from = muestra, values_from = n, values_fill = 0)
  
  mat <- as.matrix(Pfam_matrix[,-c(1,2)])  
  rownames(mat) <- Pfam_matrix$Pfam.name
  
  ##  Forzar orden alfabético
  orden_muestras <- sort(colnames(mat))
  mat <- mat[, orden_muestras]
  
  ## Chequear orden
  
  print(orden_muestras)
  
  pheatmap(mat, 
           scale = "row", 
           main = "Heatmap dominios Pfam abundantes > 10 copias (VIBRANT contigs)",
           show_rownames = TRUE,
           show_colnames = TRUE,
           cluster_cols = FALSE)
}

#Para graficar todas las proteinas encontradas
Pfam_scan

if (!is.null(Pfam_scan)) {
  Pfam_clean <- Pfam_scan %>%
    filter(!is.na(`Pfam.name`), `Pfam.name` != "") %>%
    count(muestra, `Pfam.name`, sort = TRUE) %>%
    group_by(muestra) %>%
    slice_max(n, n = 10)
  
  Pfam_matrix <- Pfam_clean %>%
    group_by(Pfam.name) %>%
    summarise(total = sum(n)) %>%
    filter(total >= 1) %>%      
    inner_join(Pfam_clean, by = "Pfam.name") %>%
    pivot_wider(names_from = muestra, values_from = n, values_fill = 0)
  
  mat <- as.matrix(Pfam_matrix[,-c(1,2)])  
  rownames(mat) <- Pfam_matrix$Pfam.name
  
  orden_muestras <- sort(colnames(mat))
  mat <- mat[, orden_muestras]
  
  ## Chequear orden
  print(orden_muestras)
  
  pheatmap(mat, 
           scale = "row", 
           main = "Heatmap dominios Pfam abundantes (VIBRANT contigs)",
           show_rownames = FALSE,
           show_colnames = TRUE,
           cluster_cols = FALSE)  # 🔑 Evita que pheatmap vuelva a reordenar
}

#Grafico 8 KO.name

Ko_name

if (!is.null(Ko_name)) {
  Ko_clean <- Ko_name %>%
    filter(!is.na(`KO.name`), `KO.name` != "") %>%
    count(muestra, `KO.name`, sort = TRUE) %>%
    group_by(muestra) %>%
    slice_max(n, n = 10)
  
  
  KO_matrix <- Ko_clean %>%
    group_by(KO.name) %>%
    summarise(total = sum(n)) %>%
    filter(total >= 10) %>%      # filtra por abundancia
    inner_join(Ko_clean, by = "KO.name") %>%
    pivot_wider(names_from = muestra, values_from = n, values_fill = 0)
  
  mat <- as.matrix(KO_matrix[,-c(1,2)])  # elimina Pfam.name y total
  rownames(mat) <- KO_matrix$KO.name
  
  orden_muestras <- sort(colnames(mat))
  mat <- mat[, orden_muestras]
  
  ## Chequear orden
  print(orden_muestras)
  
  pheatmap(mat, 
           scale = "row", 
           main = "Heatmap dominios Ko abundantes >10 copias (VIBRANT contigs)",
           show_rownames = TRUE,
           show_colnames = TRUE,
           cluster_cols = FALSE)  #  Evita que pheatmap vuelva a reordenar
}


# Ahora sí, los nodos tendrán el nombre correcto
# ------------------------------
# Resumen estadístico en consola
# ------------------------------
cat("=== RESUMEN ESTADÍSTICO ===\n")
if (!is.null(prediccion_df)) {
  cat("Predicciones:\n")
  print(table(prediccion_df$muestra, prediccion_df$prediction))
  cat("\n")
}
if (!is.null(tipo_virus_df)) {
  cat("Tipos virales:\n")
  print(table(tipo_virus_df$muestra, tipo_virus_df$type))
  cat("\n")
}
if (!is.null(vog_name_df)) {
  cat("Ejemplo de VOG names:\n")
  print(head(vog_name_df, 10))
  cat("\n")
}
if (!is.null(calidad_df)) {
  cat("Calidad de predicciones:\n")
  print(table(calidad_df$muestra, calidad_df$Quality))
  cat("\n")
}
if (!is.null(AMG_df)) {
  cat("AMGs detectados:\n")
  print(table(AMG_df$muestra, AMG_df$AMG))
  cat("\n")
}
if (!is.null(Pfam_scan)) {
  cat("Ejemplo de Pfam names:\n")
  print(head(Pfam_scan, 10))
  cat("\n")
}

"Hacer analisis con numero de KO y base de datos recien descargada de KO para poder generar el perfil metabolico con vias keeg, 
ver si se puede usar programa ssearch, y usar las vias kegg en vez de la utilizacion de las vias dadas por KO, que no son muy informativas
"

library(readr)
library(tibble)

ko_list <- read_delim("/home/jesus/CH_22/metagenoma/ko_list.txt",
                      delim = "\t", col_names = TRUE)
names(ko_list) 
names(ko_numero)
# Debe tener columnas: KO y Pathway (función/metabolismo)
# --- 2. Emparejar tus KOs con funciones ---
ko_annot <- ko_numero %>%
  left_join(ko_list %>% select(knum, definition), by = c("KO" = "knum"))

muestra_definition <- ko_annot %>%
  select(muestra, definition) %>%
  distinct()  

# Ver las primeras filas
head(muestra_definition)




library(KEGGREST)
library(pathview)
listDatabases()
library(dplyr)
library(tidyr)
library(RColorBrewer)

BiocManager::install("pathview")

colnames(ko_numero)
head(ko_numero)

# --- 1. Contar KOs por muestra ---
ko_counts <- ko_numero %>%
  group_by(muestra, KO) %>%
  summarise(Abundance = n(), .groups = "drop")

head(ko_counts)

# --- 2. Obtener todos los enlaces KO → pathway ---
all_links <- keggLink("pathway", "ko")
all_links_df <- data.frame(
  KO = sub("ko:", "", names(all_links)),
  Pathway = sub("path:", "", all_links),
  stringsAsFactors = FALSE
)

# --- 3. Filtrar solo los KO presentes en tu dataset ---
ko_to_pathway <- subset(all_links_df, KO %in% ko_counts$KO)

# --- 4. Unir abundancias con pathways ---
ko_table <- merge(ko_counts, ko_to_pathway, by = "KO")

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

# --- 2. Gráfico de barras ---
ggplot(top_pathways, aes(x = muestra, y = Abundance, fill = Description)) +
  geom_col() +
  labs(
    title = "Vía metabólica más abundante por muestra",
    x = "Muestra",
    y = "Abundancia de KO"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


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
keywords <- c("iron", "sulfate", "sulfur", "methanogenesis",
              "nitrogen", "carbon fixation", "photosynthesis",
              "fermentation", "oxidative phosphorylation", "nitrate","ammonia",
              "methane","sulfate", "anammox", "methyl","oxidation", "reduction")

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


