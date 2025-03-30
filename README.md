# Bioinformatics Project

## Prerequisites

1. **Python**: Install Python 3.8 or higher.
2. **Python Libraries**: Install the required libraries using pip:
   ```bash
   pip install numpy pandas scikit-learn biopython matplotlib requests
   ```
3. **Docker**: Install Docker to run MEME Suite in a container. Follow the [Docker installation guide](https://docs.docker.com/get-docker/).

## Setup

1. Download this repository to the local machine.
2. Download and ensure Docker is running.
3. Download the MEME Suite Docker image:
   ```bash
   docker pull memesuite/memesuite:latest
   ```

## Running the Project

1. **Prepare Input Files**:

   - Make sure an input FASTA file (e.g., `crp0.fna`) is in the project directory.

2. **Run the Script**:

   - Execute the Python script:
     ```bash
     python main.py
     ```

3. **View Results**:
   - The results will be saved in an HTML file (`results.html`) and served locally on port 8000.
   - Open the browser and navigate to:
     ```
     http://localhost:8000/results.html
     ```
