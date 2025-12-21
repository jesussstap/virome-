#!/bin/bash
# Recursively annotate assemblies with Prokka and map coverage + hmmer using user Python script
# Autor: ChatGPT + Jesus
# Uso: bash run_recursive_annotation.sh /ruta/base script_map.py

BASE_DIR="$1"
PY_SCRIPT="$2"

if [[ -z "$BASE_DIR" || -z "$PY_SCRIPT" ]]; then
    echo "Uso: bash run_recursive_annotation.sh <base_dir> <script_python>"
    exit 1
fi

# --------------------------
# 1. Buscar FASTA recursivamente
# --------------------------
mapfile -t FASTAS < <(find "$BASE_DIR" -type f \( -iname "*.fa" -o -iname "*.fna" -o -iname "*.fasta" \))

# 2. Buscar HMMER recursivamente
mapfile -t HMMERS < <(find "$BASE_DIR" -type f \( -iname "*.tsv" -o -iname "*.out" -o -iname "*.tbl" -o -iname "*.hmmer" \))

echo "Encontrados:"
echo "  FASTA  : ${#FASTAS[@]}"
echo "  HMMER  : ${#HMMERS[@]}"
echo

# --------------------------
# Función: emparejar FASTA ↔ TSV
# --------------------------
find_hmmer_match() {
    local fasta="$1"
    local base=$(basename "$fasta")
    local root="${base%.*}"

    # buscar coincidencia exacta
    for h in "${HMMERS[@]}"; do
        local hbase=$(basename "$h")
        local hroot="${hbase%.*}"
        if [[ "$hroot" == "$root" ]]; then
            echo "$h"
            return
        fi
    done

    # fallback: primer hmmer en misma carpeta
    local dir=$(dirname "$fasta")
    local fallback=$(find "$dir" -maxdepth 1 -type f \
        \( -iname "*.tsv" -o -iname "*.out" -o -iname "*.tbl" -o -iname "*.hmmer" \) | head -n 1)

    echo "$fallback"
}

# --------------------------
# Loop principal
# --------------------------
for fasta in "${FASTAS[@]}"; do
    base=$(basename "$fasta")
    sample="${base%.*}"

    echo "==============================================="
    echo "Procesando muestra: $sample"
    echo "FASTA: $fasta"

    # buscar HMMER
    hmmer=$(find_hmmer_match "$fasta")

    if [[ -z "$hmmer" ]]; then
        echo "⚠️  No se encontró archivo HMMER para $sample — saltando"
        continue
    fi

    echo "HMMER: $hmmer"

    # OUTDIR
    OUTDIR="${BASE_DIR}/ANNOT_${sample}"
    mkdir -p "$OUTDIR"

    # --------------------------
    # 3. Ejecutar PROKKA
    # --------------------------
    echo "→ Ejecutando PROKKA..."
    prokka \
        --force \
        --outdir "$OUTDIR/prokka" \
        --prefix "$sample" \
        "$fasta"

    PROKKA_GFF="$OUTDIR/prokka/${sample}.gff"
    PROKKA_FNA="$OUTDIR/prokka/${sample}.fna"

    if [[ ! -f "$PROKKA_GFF" || ! -f "$PROKKA_FNA" ]]; then
        echo "❌ ERROR: Prokka no generó archivos para $sample"
        continue
    fi

    # --------------------------
    # 4. Ejecutar TU SCRIPT PYTHON EXACTO
    # --------------------------
    echo "→ Ejecutando script Python: $PY_SCRIPT"

    python3 "$PY_SCRIPT" \
        --gff "$PROKKA_GFF" \
        --fasta "$fasta" \
        --hmmer "$hmmer" \
        --out "$OUTDIR/${sample}_annot_cov.tsv"

    echo "✔️ Completado: $sample"
    echo
done

echo "==============================================="
echo "             PIPELINE COMPLETO"
echo "==============================================="

