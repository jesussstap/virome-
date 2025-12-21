#!/bin/bash

input="/home/jesus/CH_22/galaxy_downloads/vibrant_jesus.tapia"
lecturas="/home/jesus/Lecturas_concatenadas"
salida="/home/jesus/smb_elilab/E3YTIV~H/genomics/CH_22/coverage_contigs"
threads=8

mkdir -p "$salida"

set -e
set -x  # Muestra cada comando ejecutado

# Función para procesar un contig
process_contig() {
    contig="$1"
    base=$(basename $(dirname "$contig"))


    # Generar paths de las lecturas
    R1="$lecturas/${base}_R1.fastq.gz"
    R2="$lecturas/${base}_R2.fastq.gz"
    
    sorted_bam="$salida/${base}_sorted.bam"
    coverage_out="$salida/${base}_coverage.txt"

    if [[ -f "$R1" && -f "$R2" ]]; then
        echo ">>> Procesando $base"
        
        # Mapeo + conversión + ordenación en un solo pipe
        minimap2 -t $threads -ax sr "$contig" "$R1" "$R2" \
            | samtools view -bS - \
            | samtools sort -@ $threads -o "$sorted_bam"
        
        samtools index "$sorted_bam"
        samtools coverage "$sorted_bam" > "$coverage_out"
        echo "✅ Cobertura guardada en $coverage_out"
    else
        echo "⚠️ No se encontraron las lecturas para $base"
    fi
}

export -f process_contig
export salida threads lecturas

find "$input" -type f -name "*.fasta" | while read contig; do
    process_contig "$contig"
done


