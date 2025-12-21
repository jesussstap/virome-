#!/bin/bash

# Carpeta de resultados
salida="/home/jesus/CH_22/resultados_kofamscan"
mkdir -p "$salida"

# Carpeta de entrada con subcarpetas
input="/home/jesus/CH_22/galaxy_downloads/prokka_jesus.tapia"

# Archivos de KO
ko_list="/home/jesus/CH_22/metagenoma/ko_list.txt"

# Carpeta de perfiles HMM 
profiles="/home/jesus/CH_22/metagenoma/ko_fam_hmm/ko_hmms.hmm"

# Recorrer todas las subcarpetas y archivos .fasta o .fasta.gz
find "$input" -type f \( -name "*_faa" -o -name "*_faa.fasta" \) | while read archivo; do

    # Obtener nombre del archivo sin ruta
    nombre=$(basename "$archivo")

    # Archivo de salida	
    archivo_salida="$salida/${nombre}.tsv"

    echo "🔹 Procesando $archivo ..."
    
    # Ejecutar kofamscan
    exec_annotation \
        -o "$archivo_salida" \
        -p "$profiles" \
        -k "$ko_list" \
        --cpu 2 \
        -f detail-tsv \
        "$archivo"
    
    echo "✅ Resultados guardados en $archivo_salida"

done

echo "🚀 Todos los archivos han sido procesados."

