import os
import re
import requests
import numpy as np
from Bio import pairwise2

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
        matrix = [row for row in matrix if re.match(r"^\s*\d", row)]
        sequence = ''.join(['ACGT'[np.argmax(list(map(float, row.split())))] for row in matrix])
        jaspar_motifs[motif_id] = sequence
    return jaspar_motifs

def convert_iupac_to_standard(sequence):
    iupac_to_standard = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': 'A', 'Y': 'C', 'S': 'G', 'W': 'A',
        'K': 'G', 'M': 'A', 'B': 'C', 'D': 'A',
        'H': 'A', 'V': 'A', 'N': 'A'
    }
    return ''.join([iupac_to_standard.get(base, base) for base in sequence])

def normalize_motif(motif, background_frequencies):
    normalized_motif = []
    for base in motif:
        if base in background_frequencies:
            normalized_motif.append(base)
        else:
            normalized_motif.append('A') 
    return ''.join(normalized_motif)

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
        if best_match and best_score > 0.7:
            jaspar_matches.append((meme_motif, best_match, best_score))
    return jaspar_matches
