# Bioinformatics Project

This project analyzes DNA sequences to identify motifs using MEME Suite, performs PCA and clustering, and compares motifs with the JASPAR database.

## Prerequisites

1. **Python**: Install Python 3.8 or higher.
2. **Python Libraries**: Install the required libraries using pip:
   ```bash
   pip install numpy pandas scikit-learn biopython matplotlib requests
   ```
3. **Docker**: Install Docker to run MEME Suite in a container. Follow the [Docker installation guide](https://docs.docker.com/get-docker/).

## Setup

1. Clone or download this repository to your local machine.
2. Ensure Docker is running on your system.

## Running the Project

1. **Prepare Input Files**:

   - Place your input FASTA file (e.g., `crp0.fna`) in the project directory.

2. **Run the Script**:

   - Execute the Python script:
     ```bash
     python project.py
     ```

3. **Docker Integration**:

   - The script uses Docker to run MEME Suite. Ensure Docker is running and accessible.
   - MEME Suite will be executed in a container with the following command:
     ```bash
     docker run --rm -v <project_directory>:/data memesuite/memesuite:latest
     ```
     Replace `<project_directory>` with the absolute path to your project directory.

4. **View Results**:
   - The results will be saved in an HTML file (`results.html`) and served locally on port 8000.
   - Open your browser and navigate to:
     ```
     http://localhost:8000/results.html
     ```

## Notes

- Ensure your input FASTA file is properly formatted.
- The script will automatically calculate background frequencies and fetch motifs from the JASPAR database.
- Adjust parameters like motif length and clustering range in the `main()` function of `project.py` as needed.

## Troubleshooting

- If Docker commands fail, ensure Docker is installed and running.
- If Python dependencies are missing, re-run the pip install command.
- For issues with MEME Suite, verify the input file and parameters.
