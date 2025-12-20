#ejecutar en entrono con vamb instalado 
import subprocess
import os
import glob
import shutil

def run_command(cmd):
    print(f"\nEjecutando: {cmd}\n")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise Exception(f"Fallo el comando: {cmd}")
        
def main():
    # Rutas a scripts
    concatenate_script = "~/vamb/src/concatenate.py"
    merge_aemb_script = "~/vamb/src/merge_aemb.py"
    create_fasta_script = "/vamb/src/create_fasta.py"

    # Directorios
    contig_dir = "/home/jesus/CH_22/metagenoma/ensamble_prueba"
    read_dir = "/home/jesus/Lecturas_concatenadas"
    output_dir = "/home/jesus/resultados_vamb_binsG"

    os.makedirs(output_dir, exist_ok=True)

    # Obtener archivos fasta en contig_dir
    contig_files = sorted(glob.glob(os.path.join(contig_dir, "*.fasta")))
    if len(contig_files) == 0:
        raise Exception(f"No se encontraron archivos fasta en {contig_dir}")
    try:
        min_bin_size = int(input("Ingrese el tamaño mínimo de bin en pb (ej: 6000): "))
    except ValueError:
        raise Exception("Debes ingresar un número entero válido para min_bin_size.")

    for contig in contig_files:
        sample_name = os.path.basename(contig).replace(".fasta", "")
        print(f"\nProcesando muestra: {sample_name}")

        # Archivos R1 y R2 para esta muestra
        R1 = os.path.join(read_dir, f"{sample_name}_R1.fastq.gz")
        R2 = os.path.join(read_dir, f"{sample_name}_R2.fastq.gz")

        if not os.path.exists(R1) or not os.path.exists(R2):
            print(f"⚠️ Archivos R1/R2 no encontrados para {sample_name}, saltando...")
            continue

        # Crear carpeta de output para la muestra
        sample_outdir = os.path.join(output_dir, sample_name)
        os.makedirs(sample_outdir, exist_ok=True)

        # Archivos de salida
        contigs_output = os.path.join(sample_outdir, f"{sample_name}_contigs.fna.gz")
        aemb_output = os.path.join(sample_outdir, f"{sample_name}_aemb.tsv")
        aemb_dir = os.path.join(sample_outdir, "aemb_dir")
        abundance_tsv = os.path.join(sample_outdir, f"{sample_name}_abundance.tsv")
        vamb_outdir = os.path.join(sample_outdir, "vambout")
        vae_clusters = os.path.join(vamb_outdir, "vae_clusters_split.tsv")

        os.makedirs(aemb_dir, exist_ok=True)

        # BORRAR vambout si existe
        if os.path.exists(vamb_outdir):
            print(f"Eliminando directorio existente {vamb_outdir}")
            shutil.rmtree(vamb_outdir)

        # Paso 1: Concatenar contigs (aun siendo uno, requerido por vamb)
        run_command(f"python {concatenate_script} {contigs_output} {contig}")

        # Paso 2: strobealign
        run_command(f"strobealign -t 8 --aemb {contigs_output} {R1} {R2} > {aemb_output}")

        # Paso 3: mover archivo aemb a carpeta
        run_command(f"mv {aemb_output} {aemb_dir}/")

        # Paso 4: merge aemb para tabla de abundancia
        run_command(f"python {merge_aemb_script} {aemb_dir} {abundance_tsv}")

        # Paso 5: binning con vamb
        run_command(f"vamb bin default --outdir {vamb_outdir} --fasta {contigs_output} --abundance_tsv {abundance_tsv}")

        # Paso 6: crear fasta bins
        run_command(f"python {create_fasta_script} {contigs_output} {vae_clusters} {min_bin_size} {vamb_outdir}")

    print("\n🚀  Listo.")

if __name__ == "__main__":
    main()
