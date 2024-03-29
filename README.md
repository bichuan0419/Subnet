# Subnet: Covariate-Related Functional Brain Subnetwork Extraction Tool 

## Overview
This application is designed for the extraction of functional brain subnetworks that are related to covariates of interest. It provides a user-friendly interface to simulate or import adjacency matrices representing brain connectivity and applies statistical tests to identify subnetworks.
![Screenshot 2024-01-30 134201](https://github.com/bichuan0419/Subnet/assets/43563121/ecccf3a6-9132-4eaf-9d4a-62acb3f1ec71)
![Screenshot 2024-01-30 134244](https://github.com/bichuan0419/Subnet/assets/43563121/c99febe5-acc5-480c-aae9-4fa9fa1e8918)

## Features
- **Simulate Adjacency Matrix**: Users can generate adjacency matrices based on predefined parameters such as graph structure and covariates.
- **Import Adjacency Matrix**: Users have the option to import existing adjacency matrices from `.txt` or `.csv` files.
- **Reordered Adjacency Matrix Visualization**: The application offers a visual representation of the reordered adjacency matrix to highlight the structure within the network.
- **Statistical Test**: A built-in functionality to run statistical tests using a specified alpha level and algorithm.
- **Export Options**: Users can export the results and figures generated by the application for further analysis or reporting.

## How to Use

### Simulating an Adjacency Matrix
1. Select the `Simulate Adjacency Matrix` tab.
2. Define the graph structure by specifying the number of nodes `N` and the covariates `|Vc|`.
3. Set the connection probability parameters `p0` and `p1`.
4. Determine the sample size and effect size for the simulation.
5. Click `Generate` to simulate the adjacency matrix.

### Importing an Adjacency Matrix
1. Go to the `Import Adjacency Matrix` tab.
2. Click `Import` to select and load an adjacency matrix file.

### Running Statistical Tests
1. Specify the alpha level (`α`) for the statistical test.
2. Choose the algorithm from the dropdown menu (e.g., `greedy`).
3. Enter the number of permutations `M` for the test.
4. Click `Run` to perform the statistical analysis.

## Results Interpretation
- The application will display the number of identified subnetworks along with their sizes and p-values.
- Visualizations of the observed and reordered adjacency matrices provide insight into the network's modular structure.

## Exporting Data
- Results can be exported by clicking the `Export` button, which will save the subnetworks' details.
- To export the figures, click the `Export Figure` button.

## Requirements
- To run this application, you will need [list any software, libraries, or dependencies that need to be installed].

## Installation
- Provide instructions on how to install and run the application.

## License
- State the license under which this application is released.


References:

1. Wu Q, Huang X, Culbreth AJ, Waltz JA, Hong LE, Chen S. Extracting brain disease-related connectome subgraphs by adaptive dense subgraph discovery. Biometrics. 2022 Dec;78(4):1566-1578. doi: 10.1111/biom.13537. Epub 2021 Aug 22. PMID: 34374075; PMCID: PMC10396394.
2. Shuo Chen, Yuan Zhang, Qiong Wu, Chuan Bi, Peter Kochunov, L Elliot Hong, Identifying covariate-related subnetworks for whole-brain connectome analysis, Biostatistics, 2023;, kxad007, https://doi.org/10.1093/biostatistics/kxad007
