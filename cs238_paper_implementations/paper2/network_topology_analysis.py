import pandas as pd
import networkx as nx
import os
from tqdm import tqdm

# Load the STRING interaction data (you should download this file from STRING website)
# Example link: https://stringdb-static.org/download/protein.links.v12.0/paper2/9606.protein.links.v12.0.txt.gz
# Note: Make sure to use the appropriate file path
file_path = "9606.protein.links.v12.0.txt.gz"
string_data = pd.read_csv(file_path, sep=" ", usecols=["protein1", "protein2", "combined_score"])

# Filter out interactions with low confidence scores (e.g., keep scores > 700)
high_confidence_data = string_data[string_data["combined_score"] > 700]

# Create the PPIN graph
G = nx.Graph()
print("Adding edges to the graph:")
for _, row in tqdm(high_confidence_data.iterrows(), total=high_confidence_data.shape[0]):
    G.add_edge(row["protein1"], row["protein2"])

# Calculate degree connectivity
degree_dict = dict(G.degree(G.nodes()))
print("\nDegree Connectivity:")
for gene, degree in degree_dict.items():
    print(f"{gene}: {degree}")

# Calculate betweenness centrality with progress bar
def betweenness_centrality_with_progress(G):
    betweenness = dict.fromkeys(G, 0.0)  # b[v]=0 for v in G
    nodes = list(G)
    for s in tqdm(nodes, desc="Calculating Betweenness Centrality"):
        # single source shortest paths
        S, P, sigma = _single_source_shortest_path_basic(G, s)
        # accumulation
        betweenness = _accumulate_basic(betweenness, S, P, sigma, s)
    # rescaling
    betweenness = _rescale(betweenness, len(G))
    return betweenness

def _single_source_shortest_path_basic(G, s):
    S = []
    P = {}
    for v in G:
        P[v] = []
    sigma = dict.fromkeys(G, 0.0)    # sigma[v]=0 for v in G
    sigma[s] = 1.0
    D = {}
    D[s] = 0
    Q = [s]
    while Q:   # use BFS to find shortest paths
        v = Q.pop(0)
        S.append(v)
        Dv = D[v]
        sigmav = sigma[v]
        for w in G[v]:
            if w not in D:
                Q.append(w)
                D[w] = Dv + 1
            if D[w] == Dv + 1:   # this is a shortest path, count paths
                sigma[w] += sigmav
                P[w].append(v)  # predecessors
    return S, P, sigma

def _accumulate_basic(betweenness, S, P, sigma, s):
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            betweenness[w] += delta[w]
    return betweenness

def _rescale(betweenness, n):
    if n <= 2:
        scale = None  # no normalization
    else:
        scale = 1 / ((n - 1) * (n - 2))
    if scale is not None:
        for v in betweenness:
            betweenness[v] *= scale
    return betweenness

betweenness_dict = betweenness_centrality_with_progress(G)

# Combine the results into a DataFrame
results = pd.DataFrame({
    'Gene': degree_dict.keys(),
    'Degree': degree_dict.values(),
    'Betweenness': betweenness_dict.values()
})

# Define hub genes based on top 20% degree and betweenness
degree_threshold = results['Degree'].quantile(0.8)
betweenness_threshold = results['Betweenness'].quantile(0.8)

hub_genes = results[(results['Degree'] >= degree_threshold) | (results['Betweenness'] >= betweenness_threshold)]
print("\nIdentified Hub Genes:")
print(hub_genes)

# Save the results to a CSV file
output_file_path = os.path.join(os.path.dirname(file_path), "hub_genes_output.csv")
hub_genes.to_csv(output_file_path, index=False)
print(f"\nHub genes saved to {output_file_path}")
