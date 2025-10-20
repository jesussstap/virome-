#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene map lineal de un bin con nombres de proteína
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

# -------------------------
# Leer GFF y concatenar contigs
# -------------------------
gff_file = "/home/jesus/galaxy_downloads/gff/binning_GS_F/Prokka_on_data_2:_gff.gff3"

def read_gff_linear(gff_path):
    """
    Lee un GFF3 y concatena los contigs en una sola línea.
    Devuelve DataFrame con start, end, strand y product (nombre proteína).
    """
    data = []
    contig_offsets = {}
    offset = 0

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "CDS":
                continue

            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            info = parts[8]

            # Extraer nombre de la proteína
            protein_name = "unknown"
            for field in info.split(";"):
                if field.startswith("product="):
                    protein_name = field.replace("product=", "")
                    break

            # Calcular desplazamiento si es un contig nuevo
            if contig not in contig_offsets:
                contig_offsets[contig] = offset
                offset += end  # suma longitud del contig

            start_lin = start + contig_offsets[contig]
            end_lin = end + contig_offsets[contig]

            data.append([start_lin, end_lin, strand, protein_name])

    df = pd.DataFrame(data, columns=["start", "end", "strand", "product"])
    return df

df_genes = read_gff_linear(gff_file)

# -------------------------
# Graficar bin lineal con nombre de proteína
# -------------------------
def plot_gene_map_linear(df, title="Gene map lineal"):
    fig, ax = plt.subplots(figsize=(15, 2))

    for idx, row in df.iterrows():
        start, end = row["start"], row["end"]
        color = "skyblue" if row["strand"]=="+" else "salmon"
        arrow = FancyArrow(start, 0, end-start, 0,
                           width=0.3, length_includes_head=True,
                           head_width=0.5, head_length=max((end-start)*0.05, 50),
                           color=color)
        ax.add_patch(arrow)
        # Añadir nombre de la proteína encima de la flecha
        ax.text((start+end)/2, 0.4, row["product"], rotation=90,
                ha="center", va="bottom", fontsize=6)

    ax.set_ylim(-1, 1)
    ax.set_xlim(0, df["end"].max()*1.05)
    ax.set_yticks([])
    ax.set_xlabel("Posición (bp)")
    ax.set_title(title)
    ax.grid(axis="x", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()

plot_gene_map_linear(df_genes, title="Gene map lineal - Bin 1")
