import os
import subprocess
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from Bio import SeqIO
from jaspar_utils import compare_with_jaspar 
from html_utils import save_results_to_html  
import matplotlib.pyplot as plt
import re
import requests  
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import http.server
import socketserver
import base64
from io import BytesIO

# Function to read input file and preprocess sequences
def preprocess_input(input_file):
    sequences = []
    for record in SeqIO.parse(input_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

# Function to calculate background model based on nucleotide frequencies
def calculate_background_model(input_file, background_model_file):
    nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total_count = 0

    for record in SeqIO.parse(input_file, "fasta"):
        for nucleotide in str(record.seq).upper():  
            if nucleotide in nucleotide_counts:
                nucleotide_counts[nucleotide] += 1
                total_count += 1

    if total_count == 0:
        raise ValueError("Input file contains no valid nucleotide sequences.")

    with open(background_model_file, 'w') as f:
        for nucleotide, count in nucleotide_counts.items():
            frequency = count / total_count
            f.write(f"{nucleotide} {frequency:.3f}\n")
            
# Function to run MEME using Docker
# Possible to modify the parameters for MEME
def run_meme(input_file, output_dir, motif_length, background_model):
    meme_command = [
        "docker", "run", "--rm",
        "-v", f"{os.getcwd()}:/data",
        "memesuite/memesuite:latest",
        "meme", f"/data/{input_file}",
        "-oc", f"/data/{output_dir}",
        "-mod", "zoops",  
        "-nmotifs", "10",  
        "-minw", str(motif_length[0]),
        "-maxw", str(motif_length[1]),
        "-dna",  
        "-bfile", f"/data/{background_model}",
        "-revcomp" 
    ]
    subprocess.run(meme_command)

# Helper function to reshape data for PCA and clustering 
def reshape_data(data):
    if len(data) == 0:
        return np.array(data).reshape(0, 1)
    if isinstance(data[0], list):
        return np.array(data)
    else:
        return np.array(data).reshape(-1, 1)

# Function to perform PCA and clustering on MEME output
# Also compares the motifs with JASPAR database and saves the results to an HTML file
def postprocess_output(meme_output_dir, k_range):
    motif_file = os.path.join(meme_output_dir, "meme.txt")
    motifs = []
    motif_pattern = re.compile(r"^MOTIF\s+(\S+)\s+MEME-\d+")

    # Read motifs from MEME.txt file
    with open(motif_file) as f:
        for line in f:
            match = motif_pattern.match(line)
            if match:
                motif_sequence = match.group(1)
                motifs.append(motif_sequence) 
    
    if len(motifs) < 2:
        raise ValueError("PCA requires at least two samples.")
    
    # Pad motifs to the same length
    max_length = max(len(motif) for motif in motifs)
    padded_motifs = [motif.ljust(max_length, '-') for motif in motifs]
    
    motif_matrix = np.array([[ord(char) for char in motif] for motif in padded_motifs])
    
    if motif_matrix.shape[0] < 2:
        raise ValueError("PCA requires at least two samples.")
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(motif_matrix)
    
    # Perform k-means clustering with penalty for too many clusters
    best_k = k_range[0]
    best_score = float('inf')
    for k in range(k_range[0], min(k_range[1] + 1, len(motifs))):
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(pca_result)
        score = kmeans.inertia_ + 10 * k  
        if score < best_score:
            best_score = score
            best_k = k
    
    kmeans = KMeans(n_clusters=best_k, random_state=42)
    kmeans.fit(pca_result)
    labels = kmeans.labels_

    # Compare with JASPAR and collect matches
    jaspar_matches = compare_with_jaspar(meme_output_dir)

    # Extract e-values for motifs
    e_value_pattern = re.compile(r"E-value\s*=\s*([\d\.eE+-]+)")
    e_values = []
    with open(motif_file) as f:
        for line in f:
            match = e_value_pattern.search(line)
            if match:
                e_values.append(float(match.group(1)))

    # Save results to HTML
    output_html = os.path.join("html", "results.html")  # Use os.path.join for a valid file path
    save_results_to_html(pca_result, labels, motifs, jaspar_matches, e_values, output_html)

    # Serve the HTML file on port 8000
    print(f"Serving results on http://localhost:8000/{output_html}")
    handler = http.server.SimpleHTTPRequestHandler
    with socketserver.TCPServer(("", 8000), handler) as httpd:
        httpd.serve_forever()