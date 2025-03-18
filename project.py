import os
import subprocess
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from Bio import SeqIO
from Bio.motifs import jaspar 
from Bio.motifs.jaspar.db import JASPAR5
import matplotlib.pyplot as plt
import re

# Function to run MEME using Docker
def run_meme(input_file, output_dir, motif_length, background_model):
    meme_command = [
        "docker", "run", "--rm",
        "-v", f"{os.getcwd()}:/data",
        "memesuite/memesuite:latest",
        "meme", f"/data/{input_file}",
        "-oc", f"/data/{output_dir}",
        "-mod", "zoops",  # Use ZOOPS model
        "-nmotifs", "10",  # Ensure sufficient number of motifs
        "-minw", str(motif_length[0]),
        "-maxw", str(motif_length[1]),
        "-dna",  # Specify DNA alphabet
        "-bfile", f"/data/{background_model}",
        "-revcomp"  # Allow motifs on either strand
        #"-evt", str(e_value_threshold)
    ]
    subprocess.run(meme_command)

# Function to preprocess input
def preprocess_input(input_file):
    sequences = []
    for record in SeqIO.parse(input_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

# Helper function to reshape data
def reshape_data(data):
    if len(data) == 0:
        return np.array(data).reshape(0, 1)
    if isinstance(data[0], list):
        return np.array(data)
    else:
        return np.array(data).reshape(-1, 1)

# Function to post-process output using k-means clustering
def postprocess_output(meme_output_dir, k_range):
    motif_file = os.path.join(meme_output_dir, "meme.txt")
    motifs = []
    motif_pattern = re.compile(r"^MOTIF\s+(\S+)\s+MEME-\d+")

    with open(motif_file) as f:
        for line in f:
            match = motif_pattern.match(line)
            if match:
                # Extract the motif sequence from the line
                motif_sequence = match.group(1)
                motifs.append(motif_sequence)
    
    print(f"Motifs found: {motifs}")  # Debugging line
    
    if len(motifs) < 2:
        raise ValueError("PCA requires at least two samples.")
    
    # Pad motifs to the same length
    max_length = max(len(motif) for motif in motifs)
    padded_motifs = [motif.ljust(max_length, '-') for motif in motifs]  # Use '-' for padding
    
    # Convert motifs to numerical representation
    motif_matrix = np.array([[ord(char) for char in motif] for motif in padded_motifs])
    
    if motif_matrix.shape[0] < 2:
        raise ValueError("PCA requires at least two samples.")
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(motif_matrix)
    
    # Perform k-means clustering
    best_k = k_range[0]
    best_inertia = float('inf')
    for k in range(k_range[0], min(k_range[1] + 1, len(motifs))):
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(pca_result)
        if kmeans.inertia_ < best_inertia:
            best_inertia = kmeans.inertia_
            best_k = k
    
    kmeans = KMeans(n_clusters=best_k)
    kmeans.fit(pca_result)
    labels = kmeans.labels_
    
    # Plot PCA result
    plt.scatter(pca_result[:, 0], pca_result[:, 1], c=labels)
    plt.title("PCA of Motifs with K-means Clustering")
    plt.xlabel("PCA1")
    plt.ylabel("PCA2")
    plt.show()

# Function to calculate background model
def calculate_background_model(input_file, background_model_file):
    nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total_count = 0

    for record in SeqIO.parse(input_file, "fasta"):
        for nucleotide in str(record.seq):
            if nucleotide in nucleotide_counts:
                nucleotide_counts[nucleotide] += 1
                total_count += 1

    with open(background_model_file, 'w') as f:
        for nucleotide, count in nucleotide_counts.items():
            frequency = count / total_count
            f.write(f"{nucleotide} {frequency:.3f}\n")

# Function to compare results with JASPAR
# def compare_with_jaspar(meme_output_dir, jaspar_db):
#     motif_file = os.path.join(meme_output_dir, "meme.txt")
#     meme_motifs = []
#     with open(motif_file) as f:
#         for line in f:
#             if line.startswith("MOTIF"):
#                 meme_motifs.append(line.split()[1])
    
#     jaspar_motifs = jaspar_db.fetch_motifs()
#     for meme_motif in meme_motifs:
#         for jaspar_motif in jaspar_motifs:
#             if meme_motif == jaspar_motif.matrix_id:
#                 print(f"Match found: {meme_motif} in JASPAR")

# Main function
def main():
    input_file = "crp0.fna"
    output_dir = "meme_output"
    motif_length = (6, 50)
    background_model = "background_model.txt"
    #e_value_threshold = 1000
    k_range = (2, 10)
    # jaspar_db = JASPAR5("JASPAR2020_CORE")

    sequences = preprocess_input(input_file)
    calculate_background_model(input_file, background_model)
    run_meme(input_file, output_dir, motif_length, background_model)
    postprocess_output(output_dir, k_range)
    # compare_with_jaspar(output_dir, jaspar_db)

if __name__ == "__main__":
    main()