#!/bin/bash
#para concatenar lecturas de 4 canales 
# Carpeta de entrada
LECTURAS_DIR="$HOME/Descargas/lectura"
SALIDA_DIR="$HOME/Lecturas_concatenada"

# Comprobar que existe la carpeta de entrada
if [ ! -d "$LECTURAS_DIR" ]; then
    echo "❌ No se encuentra la carpeta: $LECTURAS_DIR"
    exit 1
fi

# Crear carpeta de salida
mkdir -p "$SALIDA_DIR"

# Ir a la carpeta de lecturas
cd "$LECTURAS_DIR"

# Procesar cada archivo R1
for file in *_L001_R1_001.fastq.gz; do
    original_name="${file%%_L001_R1_001.fastq.gz}"
    prefix="${original_name:0:2}"
    char9="${original_name:8:1}"
    sample_name="${prefix}_${char9}"

    echo "📂 Procesando muestra: $original_name → $sample_name"

    # Concatenar R1
    cat ${original_name}_L00{1..4}_R1_001.fastq.gz > "$SALIDA_DIR/${sample_name}_R1.fastq.gz"
    # Concatenar R2
    cat ${original_name}_L00{1..4}_R2_001.fastq.gz > "$SALIDA_DIR/${sample_name}_R2.fastq.gz"

    echo "✅ Generados: ${sample_name}_R1.fastq.gz y ${sample_name}_R2.fastq.gz"
done

echo "🎉 Listo."

