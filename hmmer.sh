#!/bin/bash

# Ruta a la base de datos HMM
HMM_DB="/home/jesus/CH_22/metagenoma/ko_fam_hmm/all_hmms_clean.hmm"

# Carpeta raíz con subcarpetas de Prokka
FASTA_DIR="/home/jesus/CH_22/galaxy_downloads/PROKKA_Contigs_vir"

# Carpeta de salida
OUT_DIR="smb://elilab.local/public/E3YTIV~H/genomics/CH_22"
mkdir -p "$OUT_DIR"

# Contador total de archivos
total=$(find "$FASTA_DIR" -type f -name "*_faa.fasta" | wc -l)
i=1

# Iterar sobre todos los archivos .faa dentro de las subcarpetas
find "$FASTA_DIR" -type f -name "*_faa.fasta" | while read -r SEQ_FILE; do
    FILENAME=$(basename "$SEQ_FILE" .faa)
    PARENT_DIR=$(basename "$(dirname "$SEQ_FILE")")
    OUTPUT_SUBDIR="$OUT_DIR/$PARENT_DIR"
    mkdir -p "$OUTPUT_SUBDIR"

    DOMTBL_OUT="$OUTPUT_SUBDIR/${FILENAME}_hmmer.tabular"
    FULL_OUT="$OUTPUT_SUBDIR/${FILENAME}_hmmer.txt"

    echo "[$i/$total] Procesando: $PARENT_DIR/$FILENAME ..."
    # Verificar que el archivo tenga secuencias FASTA válidas
    if grep -q "^>" "$SEQ_FILE"; then
        hmmscan --cpu 2 \
            --domtblout "$DOMTBL_OUT" \
            "$HMM_DB" "$SEQ_FILE" > "$FULL_OUT"
        echo "  HMMER completado para $FILENAME"
    else
        echo " Archivo sin formato FASTA: $SEQ_FILE"
    fi

    i=$((i+1))
done

echo "Análisis HMMER finalizado. Resultados en: $OUT_DIR"

