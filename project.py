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
import http.server
import socketserver
import base64
from io import BytesIO

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

# Function to save plot and cluster information to an HTML file
def save_results_to_html(pca_result, labels, motifs, jaspar_matches, output_html):
    # Save the scatter plot as a base64-encoded image
    plt.figure()
    scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], c=labels, cmap="viridis")
    plt.title("PCA of Motifs with K-means Clustering")
    plt.xlabel("PCA1")
    plt.ylabel("PCA2")
    # Add a legend to show cluster numbers
    legend1 = plt.legend(*scatter.legend_elements(), title="Clusters")
    plt.gca().add_artist(legend1)
    buffer = BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode("utf-8")
    buffer.close()
    plt.close()

    # Read the HTML template
    template_path = os.path.join(os.getcwd(), "results_template.html")
    with open(template_path, "r") as template_file:
        html_template = template_file.read()

    # Generate the table rows for motifs and clusters
    table_rows = "".join(f"<tr><td>{motif}</td><td>{label}</td></tr>" for motif, label in zip(motifs, labels))

    # Generate the table rows for JASPAR matches
    jaspar_rows = "".join(
        f"<tr><td>{motif}</td><td>{match}</td><td>{score:.2f}</td></tr>"
        for motif, match, score in jaspar_matches
    )

    # Replace placeholders in the template
    html_content = (
        html_template.replace("{{scatter_plot}}", image_base64)
        .replace("{{table_rows}}", table_rows)
        .replace("{{jaspar_rows}}", jaspar_rows)
    )

    # Write the final HTML content to the output file
    with open(output_html, "w") as output_file:
        output_file.write(html_content)

# Modify postprocess_output to save results to HTML and serve it
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

    # Save results to HTML
    output_html = "results.html"
    save_results_to_html(pca_result, labels, motifs, jaspar_matches, output_html)

    # Serve the HTML file on port 8000
    print(f"Serving results on http://localhost:8000/{output_html}")
    handler = http.server.SimpleHTTPRequestHandler
    with socketserver.TCPServer(("", 8000), handler) as httpd:
        httpd.serve_forever()

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

# Modify compare_with_jaspar to return matches
def compare_with_jaspar(meme_output_dir):
    motif_file = os.path.join(meme_output_dir, "meme.txt")
    meme_motifs = []
    with open(motif_file) as f:
        for line in f:
            if line.startswith("MOTIF"):
                meme_motifs.append(line.split()[1])
    
    jaspar_motifs = fetch_jaspar_motifs()
    jaspar_background_frequencies = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    jaspar_matches = []
    for meme_motif in meme_motifs:
        standard_motif = convert_iupac_to_standard(meme_motif)
        normalized_motif = normalize_motif(standard_motif, jaspar_background_frequencies)
        best_match = None
        best_score = 0
        for jaspar_motif in jaspar_motifs.values():
            alignments = pairwise2.align.globalxx(normalized_motif, jaspar_motif)
            best_alignment = alignments[0]
            score = best_alignment[2]
            normalized_score = score / max(len(normalized_motif), len(jaspar_motif))
            if normalized_score > best_score:
                best_score = normalized_score
                best_match = jaspar_motif
        if best_match and best_score > 0.7:  # Threshold for similarity
            jaspar_matches.append((meme_motif, best_match, best_score))
    return jaspar_matches

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