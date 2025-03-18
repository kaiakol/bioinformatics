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
import requests  
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

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
    
    # Print PCA results and cluster labels
    for i, (x, y) in enumerate(pca_result):
        print(f"Motif: {motifs[i]}, PCA1: {x:.2f}, PCA2: {y:.2f}, Cluster: {labels[i]}")
    
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

# Function to fetch motifs from JASPAR database
def fetch_jaspar_motifs():
    url = "https://jaspar.genereg.net/api/v1/matrix/?format=meme"
    response = requests.get(url)
    response.raise_for_status()
    meme_format = response.text
    jaspar_motifs = {}
    motif_pattern = re.compile(r"MOTIF\s+(\S+)\s+.*?\n(.*?)\n\n", re.DOTALL)
    for match in motif_pattern.finditer(meme_format):
        motif_id = match.group(1)
        matrix = match.group(2).strip().split('\n')
        # Filter out lines that do not contain the letter-probability matrix
        matrix = [row for row in matrix if re.match(r"^\s*\d", row)]
        sequence = ''.join(['ACGT'[np.argmax(list(map(float, row.split())))] for row in matrix])
        jaspar_motifs[motif_id] = sequence
    return jaspar_motifs

# Function to convert IUPAC codes to standard nucleotide representations
def convert_iupac_to_standard(sequence):
    iupac_to_standard = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': 'A',  # A or G (purine)
        'Y': 'C',  # C or T (pyrimidine)
        'S': 'G',  # G or C
        'W': 'A',  # A or T
        'K': 'G',  # G or T
        'M': 'A',  # A or C
        'B': 'C',  # C, G, or T
        'D': 'A',  # A, G, or T
        'H': 'A',  # A, C, or T
        'V': 'A',  # A, C, or G
        'N': 'A'   # Any base (A, C, G, or T)
    }
    return ''.join([iupac_to_standard.get(base, base) for base in sequence])

# Function to normalize motif using background frequencies
def normalize_motif(motif, background_frequencies):
    normalized_motif = []
    for base in motif:
        if base in background_frequencies:
            normalized_motif.append(base)
        else:
            normalized_motif.append('A')  # Default to 'A' for unknown bases
    return ''.join(normalized_motif)

# Function to compare results with JASPAR using sequence alignment
def compare_with_jaspar(meme_output_dir):
    motif_file = os.path.join(meme_output_dir, "meme.txt")
    meme_motifs = []
    with open(motif_file) as f:
        for line in f:
            if line.startswith("MOTIF"):
                meme_motifs.append(line.split()[1])
    
    jaspar_motifs = fetch_jaspar_motifs()
    jaspar_background_frequencies = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    print(f"JASPAR motifs: {list(jaspar_motifs.values())}")
    found_match = False
    for meme_motif in meme_motifs:
        standard_motif = convert_iupac_to_standard(meme_motif)
        normalized_motif = normalize_motif(standard_motif, jaspar_background_frequencies)
        print(f"Comparing motif: {normalized_motif}")
        for jaspar_motif in jaspar_motifs.values():
            alignments = pairwise2.align.globalxx(normalized_motif, jaspar_motif)
            best_alignment = alignments[0]
            score = best_alignment[2]
            normalized_score = score / max(len(normalized_motif), len(jaspar_motif))
            if normalized_score > 0.7:  # Threshold for similarity
                print(f"Match found: {normalized_motif} in JASPAR with alignment score {score} (normalized score: {normalized_score:.2f})")
                print(format_alignment(*best_alignment))
                found_match = True
                break
    if not found_match:
        print("No matches found between MEME motifs and JASPAR motifs.")

# Main function
def main():
    input_file = "crp0.fna"
    output_dir = "meme_output"
    motif_length = (6, 50)
    background_model = "background_model.txt"
    #e_value_threshold = 1000
    k_range = (2, 10)

    sequences = preprocess_input(input_file)
    calculate_background_model(input_file, background_model)
    run_meme(input_file, output_dir, motif_length, background_model)
    postprocess_output(output_dir, k_range)
    compare_with_jaspar(output_dir)  

if __name__ == "__main__":
    main()