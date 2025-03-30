import os
import base64
import matplotlib.pyplot as plt
from io import BytesIO

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
    template_path = os.path.join(os.getcwd(), "html", "results_template.html")
    with open(template_path, "r") as template_file:
        html_template = template_file.read()

    # Sort motifs and clusters by cluster number
    sorted_motifs_clusters = sorted(zip(motifs, labels), key=lambda x: x[1])

    # Generate the table rows for motifs and clusters
    table_rows = "".join(f"<tr><td>{motif}</td><td>{label}</td></tr>" for motif, label in sorted_motifs_clusters)

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
