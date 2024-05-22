import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

# Load the hub genes CSV file
file_path = "hub_genes_output.csv"
hub_genes = pd.read_csv(file_path)

# Plot Degree Distribution
plt.figure(figsize=(10, 6))
plt.hist(hub_genes['Degree'], bins=30, color='skyblue', edgecolor='black')
plt.title('Degree Distribution')
plt.xlabel('Degree')
plt.ylabel('Frequency')
degree_plot_path = "degree_distribution.png"
plt.savefig(degree_plot_path)
print(f"Degree distribution plot saved as {degree_plot_path}")
plt.close()

# Plot Betweenness Centrality Distribution
plt.figure(figsize=(10, 6))
plt.hist(hub_genes['Betweenness'], bins=30, color='lightcoral', edgecolor='black')
plt.title('Betweenness Centrality Distribution')
plt.xlabel('Betweenness Centrality')
plt.ylabel('Frequency')
betweenness_plot_path = "betweenness_centrality_distribution.png"
plt.savefig(betweenness_plot_path)
print(f"Betweenness centrality distribution plot saved as {betweenness_plot_path}")
plt.close()

# Create a Network Visualization highlighting hub genes
# Load the STRING interaction data again to recreate the graph
string_file_path = "9606.protein.links.v12.0.txt.gz"
print("Loading STRING interaction data:")
string_data = pd.read_csv(string_file_path, sep=" ", usecols=["protein1", "protein2", "combined_score"])
high_confidence_data = string_data[string_data["combined_score"] > 700]

# Create the PPIN graph
G = nx.Graph()
print("Adding edges to the graph:")
for _, row in tqdm(high_confidence_data.iterrows(), total=high_confidence_data.shape[0]):
    G.add_edge(row["protein1"], row["protein2"])

# Highlight hub genes
hub_gene_list = hub_genes['Gene'].tolist()
node_colors = ['red' if node in hub_gene_list else 'blue' for node in tqdm(G.nodes(), desc="Coloring nodes")]

plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.1)  # Position nodes using the Fruchterman-Reingold force-directed algorithm
nx.draw_networkx_nodes(G, pos, node_size=20, node_color=node_colors, alpha=0.6)
nx.draw_networkx_edges(G, pos, alpha=0.2)
plt.title('Network Visualization Highlighting Hub Genes')
network_plot_path = "network_visualization.png"
plt.savefig(network_plot_path)
print(f"Network visualization plot saved as {network_plot_path}")
plt.close()
