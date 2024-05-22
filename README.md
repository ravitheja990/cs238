# CS238 Paper Implementations

This repository contains the implementation of methods from the paper "A Computational Framework for the Prioritization of Disease-Gene Candidates".

## Contents

- `impact_project.py`: The main script that processes the STRING interaction data, identifies hub genes, and generates various plots.
- `hub_genes_output.csv`: The output file containing the identified hub genes with their degree and betweenness centrality.
- `9606.protein.links.v12.0.txt.gz`: STRING interaction data for Homo sapiens (to be downloaded separately).
- `README.md`: This readme file.
- Other scripts and files related to the project.

## Setup

1. **Clone the Repository**

   ```bash
   git clone https://github.com/ravitheja990/cs238
   cd cs238_paper_implementations

    Install Dependencies

    Make sure you have Python installed. Install the required Python packages using pip:

    bash

pip install pandas networkx matplotlib tqdm

Download STRING Interaction Data

Download the STRING interaction file for Homo sapiens and place it in the repository directory.

bash

    wget https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz

Usage

    Run the Script

    Execute the main script to process the data and generate plots:

    bash

    python3 impact_project.py

    Generated Output
        hub_genes_output.csv: Contains the identified hub genes.
        degree_distribution.png: Plot showing the degree distribution of the genes.
        betweenness_centrality_distribution.png: Plot showing the betweenness centrality distribution of the genes.
        network_visualization.png: Network visualization highlighting the hub genes.

License

This project is licensed under the MIT License. See the LICENSE file for more details.
Acknowledgements

    The implementation is based on methods described in the paper "A Computational Framework for the Prioritization of Disease-Gene Candidates".
    STRING database for providing the interaction data.

Contact

For any questions or issues, please contact Ravi Theja at rtang069@ucr.edu.

markdown


Replace `https://github.com/ravitheja990/cs238` with the URL of your repository, and `Ravi Theja` and `rtang069@ucr.edu` with your actual name and email address.

### Steps to Use

1. **Copy the Content**: Copy the entire content above.
2. **Create the README.md File**: Create a new file named `README.md` in the root directory of your repository.
3. **Paste the Content**: Paste the copied content into the `README.md` file and save it.

This will provide a comprehensive and formatted README file for your Git repository. If you need any additional modifications, feel free to ask!

