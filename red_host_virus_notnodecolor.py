#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 10:47:19 2025

@author: jesus
"""

import os
import pandas as pd
import glob
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.community import greedy_modularity_communities
import seaborn as sns
# ==============================
# 1. Leer resultados de iPHoP
# ==============================
ruta_ihpop = "/home/jesus/CH_22/galaxy_downloads/iphop/todos_virus/todos_los_virus"
archivos_csv = glob.glob(os.path.join(ruta_ihpop, "*.csv"))

lista_dfs = []
for archivo in archivos_csv:
    df = pd.read_csv(archivo, comment="#", sep=None, engine="python")
    df["sample"] = os.path.splitext(os.path.basename(archivo))[0]
    lista_dfs.append(df)

hospedador = pd.concat(lista_dfs, ignore_index=True)
hospedador = hospedador[["Virus", "Host taxonomy", "Confidence score", "sample"]]

# ==============================
# 2. Elegir hospedador más confiable
# ==============================
confianza_mas_alta = (
    hospedador.sort_values("Confidence score", ascending=False)
    .groupby("Virus", as_index=False)
    .first()
)

# ==============================
# 3. Leer resultados de VIBRANT
# ==============================
ruta_vibrant = "/home/jesus/CH_22/galaxy_downloads/vibrant_jesus.tapia"
archivos_tabular = glob.glob(os.path.join(ruta_vibrant, "**/*.tabular"), recursive=True)
archivos_no_vacios = [f for f in archivos_tabular if os.path.getsize(f) > 0]

lista_vibrant = []
for f in archivos_no_vacios:
    df = pd.read_csv(f, sep="\t")
    df["muestra"] = os.path.basename(os.path.dirname(f))
    lista_vibrant.append(df)

vibrant = pd.concat(lista_vibrant, ignore_index=True)
vibrant = vibrant[["scaffold", "muestra", "type"]]
vibrant.rename(columns={"scaffold": "Virus"}, inplace=True)


### para graficar solo lisogenicos 
vibrant = vibrant[vibrant["type"].str.contains("lytic", case=False, na=False)]

# ==============================
# 4. Unir iPHoP + VIBRANT
# ==============================
df_unido = pd.merge(vibrant, confianza_mas_alta, on="Virus", how="inner")

# ==============================
# 5. Filtrar confianza >= 90
# ==============================
mas_confiables = df_unido[df_unido["Confidence score"] >= 90]

# ==============================
# 6. Extraer nivel taxonómico deseado
# ==============================
nivel_objetivo = "p__"  # Cambia si quieres otro nivel taxonómico

def extraer_nombre_host(tax):
    if pd.isna(tax):
        return None
    partes = [p.strip() for p in str(tax).split(";") if p.strip()]
    for p in partes:
        if p.startswith(nivel_objetivo):
            return p.split("__")[-1]
    return None

mas_confiables["Host taxonomy"] = mas_confiables["Host taxonomy"].apply(extraer_nombre_host)

# ==============================
# 7. Crear red
# ==============================
edges = mas_confiables[["Virus", "Host taxonomy", "muestra", "type"]].dropna()
G = nx.Graph()

for _, row in edges.iterrows():
    virus = row["Virus"]
    host = row["Host taxonomy"]
    muestra = row["muestra"]
    tipo = row["type"]
    G.add_node(virus, tipo="virus", muestra=muestra, tipo_viral=tipo)
    G.add_node(host, tipo="host")
    G.add_edge(virus, host)

# ==============================
# 8. Colores por muestra
# ==============================
muestras = sorted(list(set(nx.get_node_attributes(G, "muestra").values())))
n_muestras = len(muestras)
# Si hay pocas muestras usa tab20, si hay muchas usa husl (seaborn)
if n_muestras <= 20:
    base_cmap = plt.cm.get_cmap("tab20", n_muestras)
    colors = [base_cmap(i) for i in range(n_muestras)]
else:
    # husl genera colores equidistantes en el espacio perceptual HSL
    colors = sns.color_palette("husl", n_muestras)

# Crear diccionario muestra → color
color_mapa = {m: colors[i] for i, m in enumerate(muestras)}
# ==============================
# 9. Separar nodos
# ==============================
virus_nodes = [n for n, d in G.nodes(data=True) if d.get("tipo") == "virus"]
host_nodes = [n for n, d in G.nodes(data=True) if d.get("tipo") == "host"]

# Asignar color según muestra
virus_colors = [color_mapa[G.nodes[n]["muestra"]] for n in virus_nodes]

# ==============================
# 10. Asignar forma según tipo viral
# ==============================
tipos_virales = list(sorted(set([G.nodes[n]["tipo_viral"] for n in virus_nodes])))
shapes = ["^", "s"]  # diferentes formas disponibles
shape_map = {t: shapes[i % len(shapes)] for i, t in enumerate(tipos_virales)}

# ==============================
# 11. Layout y comunidades
# ==============================
pos = nx.spring_layout(G, k=1, iterations=300, seed=42)
comunidades = list(greedy_modularity_communities(G))
comunidad_map = {}
for i, c in enumerate(comunidades):
    for node in c:
        comunidad_map[node] = i

# ==============================
# 12. Detectar hosts compartidos entre muestras
# ==============================
host_multi = []
for host in host_nodes:
    muestras_conectadas = set(edges.loc[edges["Host taxonomy"] == host, "muestra"])
    if len(muestras_conectadas) > 1:
        host_multi.append(host)
        
# ==============================
# 12b. Crear tabla con especies presentes en más de una muestra
# ==============================
host_multi_info = (
    edges.groupby("Host taxonomy")["muestra"]
    .apply(lambda x: sorted(set(x)))
    .reset_index()
)

# Filtrar solo los hosts que aparecen en más de una muestra
host_multi_info = host_multi_info[host_multi_info["muestra"].apply(len) > 1]

print("\n📊 Especies presentes en más de una muestra:")
print(host_multi_info.to_string(index=False))

# Si quieres exportarlo a CSV
host_multi_info.to_csv("hosts_compartidos.csv", index=False)
# ==============================
# 13. Dibujar red
# ==============================
plt.figure(figsize=(22, 16))

# Aristas
nx.draw_networkx_edges(G, pos, alpha=0.5, edge_color="black")

# Hosts normales (negros)
nx.draw_networkx_nodes(G, pos,
                       nodelist=[h for h in host_nodes if h not in host_multi],
                       node_color="gray",
                       node_shape="o",
                       node_size=500,
                       alpha=0.5,
                       linewidths=0.5,
                       edgecolors="black")

# Hosts compartidos (negros con borde rojo)
nx.draw_networkx_nodes(G, pos,
                       nodelist=host_multi,
                       node_color="black",
                       node_shape="o",
                       node_size=500,
                       alpha=0.7,
                       linewidths=2.5,
                       edgecolors="red")

# Virus (color por muestra, forma por tipo)
for tipo, shape in shape_map.items():
    virus_tipo = [n for n in virus_nodes if G.nodes[n]["tipo_viral"] == tipo]
    nx.draw_networkx_nodes(G, pos,
                           nodelist=virus_tipo,
                           node_color=[color_mapa[G.nodes[n]["muestra"]] for n in virus_tipo],
                           node_shape=shape,
                           node_size=180,
                           alpha=0.85,
                           linewidths=1,
                           edgecolors="gray")
                          

# Etiquetas de hosts
label_pos = {n: (x, y + 0.05) for n, (x, y) in pos.items() if n in host_nodes}
nx.draw_networkx_labels(G, label_pos,
                        labels={n: n for n in host_nodes},
                        font_size=19,
                        font_color="black",
                        font_weight="bold")

# ==============================
# 14. Leyenda
# ==============================
# Colores por muestra
for m, c in color_mapa.items():
    plt.scatter([], [], c=[c], label=m, marker="o", s=180)
# Formas por tipo viral
for tipo, shape in shape_map.items():
    plt.scatter([], [], c="white", edgecolor="black", label=f"{tipo} virus", marker=shape, s=250, linewidth=1.5)
# Hosts compartidos
plt.scatter([], [], c="black", edgecolor="red", linewidth=2.5,
            label="Shared host", marker="o", s=250)

plt.legend(title="Samples and viral type", bbox_to_anchor=(1.05, 1),
           loc="upper left", fontsize=20, title_fontsize=25)
plt.title(f"Virus–host network (Level {nivel_objetivo}, Confidence ≥ 90%)", fontsize=40)
plt.axis("off")
plt.tight_layout()
plt.show()
