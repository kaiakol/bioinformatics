import os
import base64
import matplotlib.pyplot as plt
from io import BytesIO

# Function to save clustering and JASPAR results to HTML file
def save_results_to_html(pca_result, labels, motifs, e_values, tomtom_matches, output_html):
    plt.figure()
    scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], c=labels, cmap="viridis")
    plt.title("PCA of Motifs with K-means Clustering")
    plt.xlabel("PCA1")
    plt.ylabel("PCA2")
    legend1 = plt.legend(*scatter.legend_elements(), title="Clusters")
    plt.gca().add_artist(legend1)
    buffer = BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode("utf-8")
    buffer.close()
    plt.close()

    template_path = os.path.join(os.getcwd(), "html", "results_template.html")
    with open(template_path, "r") as template_file:
        html_template = template_file.read()

    # Prepare table rows with motif, cluster, and e-value
    sorted_motifs_clusters = sorted(zip(motifs, labels, e_values), key=lambda x: x[1])
    table_rows = "".join(
        f"<tr><td>{motif}</td><td>{label}</td><td>{e_value:.2e}</td></tr>"
        for motif, label, e_value in sorted_motifs_clusters
    )

    # Prepare table rows with TOMTOM matches
    tomtom_rows = "".join(
        f"<tr><td>{meme_motif}</td><td>{jaspar_motif}</td><td>{evalue:.2e}</td></tr>"
        for meme_motif, jaspar_motif, evalue in tomtom_matches
    )

    html_content = (
        html_template.replace("{{scatter_plot}}", image_base64)
        .replace("{{table_rows}}", table_rows)
        .replace("{{tomtom_rows}}", tomtom_rows)
    )

    with open(output_html, "w") as output_file:
        output_file.write(html_content)
