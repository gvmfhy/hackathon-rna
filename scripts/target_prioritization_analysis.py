import pandas as pd
import numpy as np
import os
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import gzip
import pyarrow.parquet as pq
from datetime import datetime

print("Starting miRNA Target Prioritization Analysis...")

# Get current directory and set paths
current_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(current_dir)
results_dir = os.path.join(base_dir, 'results')
data_dir = os.path.join(base_dir, 'data')
string_dir = os.path.join(data_dir, 'string_db')
ot_dir = os.path.join(data_dir, 'open_targets')

# Create output directories
network_dir = os.path.join(results_dir, 'networks')
os.makedirs(network_dir, exist_ok=True)

# Current timestamp for output filenames
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Set file paths
mirna_significant_file = os.path.join(results_dir, 'mirna_de_significant_with_names.csv')
# Point to the new mirnet subdirectory for miRTarBase data
mirtarbase_file = os.path.join(data_dir, 'mirnet', 'miRNet-mir-gene-hsa-mirtarbase.csv')
string_links_file = os.path.join(string_dir, '9606.protein.links.v11.5.txt.gz')
string_info_file = os.path.join(string_dir, '9606.protein.info.v11.5.txt.gz')

# Verify if files exist
for file_path in [mirna_significant_file, mirtarbase_file, string_links_file, string_info_file]:
    if not os.path.exists(file_path):
        print(f"Error: File not found - {file_path}")
        exit(1)

print("All required input files found.")

# Configuration parameters
STRING_SCORE_THRESHOLD = 700  # Confidence score threshold for STRING PPIs (700 = high confidence)
HUB_CENTRALITY_PERCENTILE = 90  # Top X percentile of genes to consider as hubs

# -----------------------------------------
# 1. Load and Process Data
# -----------------------------------------

# Load miRNA significant results file
print("Loading significantly differentially expressed miRNAs...")
df_mirna = pd.read_csv(mirna_significant_file)
print(f"Number of significant miRNAs: {len(df_mirna)}")

# Load the miRTarBase file
print("Loading miRTarBase miRNA-target interactions...")
df_targets = pd.read_csv(mirtarbase_file)
print(f"Number of miRNA-target interactions in miRTarBase: {len(df_targets)}")

# Create a standardized lowercase version for matching miRNAs
df_mirna['mirna_lower'] = df_mirna['miRNA_name'].str.lower().str.replace('mir-', 'mir-')
df_targets['mirna_lower'] = df_targets['mir_id'].str.lower().str.replace('mir-', 'mir-')

# Get lists of upregulated and downregulated miRNAs
up_mirnas = df_mirna[df_mirna['log2FC'] > 0]['mirna_lower'].tolist()
down_mirnas = df_mirna[df_mirna['log2FC'] < 0]['mirna_lower'].tolist()

print(f"Upregulated miRNAs: {up_mirnas}")
print(f"Downregulated miRNAs: {down_mirnas}")

# Filter for targets of upregulated miRNAs
up_targets_df = df_targets[df_targets['mirna_lower'].isin(up_mirnas)]
up_target_genes = up_targets_df['symbol'].unique().tolist()

# Filter for targets of downregulated miRNAs
down_targets_df = df_targets[df_targets['mirna_lower'].isin(down_mirnas)]
down_target_genes = down_targets_df['symbol'].unique().tolist()

print(f"Number of target genes for upregulated miRNAs: {len(up_target_genes)}")
print(f"Number of target genes for downregulated miRNAs: {len(down_target_genes)}")

# -----------------------------------------
# 2. Load and Process STRING Data
# -----------------------------------------
print("Loading STRING protein-protein interaction data...")

# Load STRING protein info
print("Loading STRING protein info...")
string_info = {}
with gzip.open(string_info_file, 'rt') as f:
    # Skip header
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        string_id = parts[0]  # String ID (e.g., 9606.ENSP00000000233)
        gene_name = parts[1]  # Gene name/symbol
        string_info[string_id] = gene_name

print(f"Loaded {len(string_info)} protein ID to gene symbol mappings from STRING")

# Load STRING PPI data with high confidence score
print("Loading STRING protein-protein interactions (this may take a moment)...")
interactions = []
with gzip.open(string_links_file, 'rt') as f:
    # Skip header
    next(f)
    for line in f:
        protein1, protein2, score = line.strip().split()
        score = int(score)
        if score >= STRING_SCORE_THRESHOLD:
            interactions.append((protein1, protein2, score))

print(f"Loaded {len(interactions)} high-confidence protein interactions from STRING")

# Function to convert STRING IDs to gene symbols
def string_id_to_symbol(string_id):
    return string_info.get(string_id, string_id)

# -----------------------------------------
# 3. Load Open Targets Data for Druggability Assessment
# -----------------------------------------
print("Loading Open Targets target data for druggability assessment...")

# Function to read and process parquet files
def read_open_targets_data():
    all_data = []
    # Get all parquet files in the directory
    parquet_files = [f for f in os.listdir(ot_dir) if f.endswith('.parquet')]
    
    print(f"Processing {len(parquet_files)} Open Targets parquet files...")
    for i, file in enumerate(parquet_files):
        if i % 10 == 0:  # Progress report every 10 files
            print(f"  Processing file {i+1}/{len(parquet_files)}: {file}")
        
        file_path = os.path.join(ot_dir, file)
        try:
            # Read parquet file
            table = pq.read_table(file_path)
            df = table.to_pandas()
            
            # Extract relevant columns for druggability assessment
            # We need gene ID (id), gene symbol (approvedSymbol), and tractability info
            if 'id' in df.columns and 'approvedSymbol' in df.columns and 'tractability' in df.columns:
                # Select only needed columns to save memory
                subset = df[['id', 'approvedSymbol', 'tractability']]
                all_data.append(subset)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    # Combine all dataframes
    if all_data:
        return pd.concat(all_data, ignore_index=True)
    else:
        return pd.DataFrame(columns=['id', 'approvedSymbol', 'tractability'])

# Load Open Targets data
ot_data = read_open_targets_data()
print(f"Loaded {len(ot_data)} target entries from Open Targets")

# Function to assess druggability from Open Targets data
def assess_druggability(gene_symbol):
    """
    Assess if a gene is considered druggable based on Open Targets data.
    Returns a tuple containing:
    (is_druggable, druggability_info)
    """
    # Initialize default druggability info dictionary
    druggability_info = {'small_molecule': False, 'antibody': False, 'other': False}
    
    # Find the gene in Open Targets data
    gene_records = ot_data[ot_data['approvedSymbol'] == gene_symbol]
    
    if len(gene_records) == 0:
        return (False, druggability_info)
    
    # Check tractability information
    is_druggable = False
    
    try:
        for _, record in gene_records.iterrows():
            # Check if tractability data exists and is not None/empty
            tractability_data = record.get('tractability')
            if tractability_data is not None and len(tractability_data) > 0:
                # Iterate through the items (should work for lists and numpy arrays)
                for tract_item in tractability_data:
                    # Ensure the item is a dictionary before accessing keys
                    if isinstance(tract_item, dict):
                        modality = tract_item.get('modality', '')
                        value = tract_item.get('value', False)
                        
                        if value:  # If the tractability value is true
                            is_druggable = True
                            
                            modality_lower = modality.lower()
                            if 'small' in modality_lower:
                                druggability_info['small_molecule'] = True
                            elif 'antibody' in modality_lower:
                                druggability_info['antibody'] = True
                            else:
                                # Capture any other modality marked as True
                                if modality: # Avoid counting empty modalities
                                    druggability_info['other'] = True
    except Exception as e:
        print(f"Warning: Error processing druggability for {gene_symbol}: {e}")
        # Continue with default values
    
    return (is_druggable, druggability_info)

# -----------------------------------------
# 4. Build Target Gene Networks
# -----------------------------------------
print("Building target gene networks...")

# Create networks for up and down regulated miRNA targets
G_up_targets = nx.Graph()
G_down_targets = nx.Graph()

# Add nodes (genes) to the networks
for gene in up_target_genes:
    G_up_targets.add_node(gene)

for gene in down_target_genes:
    G_down_targets.add_node(gene)

# Add edges (protein-protein interactions) to the networks
print("Adding protein-protein interactions to the networks...")
for protein1, protein2, score in interactions:
    gene1 = string_id_to_symbol(protein1)
    gene2 = string_id_to_symbol(protein2)
    
    # Normalize the score between 0 and 1
    norm_score = score / 1000
    
    # Add edge to up-regulated miRNA targets network
    if gene1 in up_target_genes and gene2 in up_target_genes:
        G_up_targets.add_edge(gene1, gene2, weight=norm_score)
    
    # Add edge to down-regulated miRNA targets network
    if gene1 in down_target_genes and gene2 in down_target_genes:
        G_down_targets.add_edge(gene1, gene2, weight=norm_score)

print(f"Up-regulated miRNA targets network: {G_up_targets.number_of_nodes()} nodes, {G_up_targets.number_of_edges()} edges")
print(f"Down-regulated miRNA targets network: {G_down_targets.number_of_nodes()} nodes, {G_down_targets.number_of_edges()} edges")

# -----------------------------------------
# 5. Network Analysis (Centrality & Hub Identification)
# -----------------------------------------
print("Analyzing networks for centrality and hub identification...")

# Calculate centrality measures and identify hub genes for up-regulated miRNA targets
def analyze_network(G, name):
    if G.number_of_nodes() == 0:
        print(f"No nodes in {name} network. Skipping analysis.")
        return pd.DataFrame()
    
    centrality_measures = {}
    
    # Calculate degree centrality
    print(f"Calculating degree centrality for {name} network...")
    centrality_measures['degree'] = nx.degree_centrality(G)
    
    # Calculate betweenness centrality (only if there are edges)
    if G.number_of_edges() > 0:
        print(f"Calculating betweenness centrality for {name} network...")
        centrality_measures['betweenness'] = nx.betweenness_centrality(G)
    else:
        centrality_measures['betweenness'] = {node: 0.0 for node in G.nodes()}
    
    # Create a DataFrame with centrality measures
    centrality_df = pd.DataFrame(index=G.nodes())
    for measure, values in centrality_measures.items():
        centrality_df[measure] = pd.Series(values)
    
    # Calculate combined centrality score (average of normalized measures)
    centrality_df['combined_score'] = centrality_df.mean(axis=1)
    
    # Sort by combined score
    centrality_df = centrality_df.sort_values('combined_score', ascending=False)
    
    # Identify hub genes (top percentile)
    hub_threshold = np.percentile(centrality_df['combined_score'].values, HUB_CENTRALITY_PERCENTILE)
    centrality_df['is_hub'] = centrality_df['combined_score'] >= hub_threshold
    
    # Add druggability information
    print(f"Assessing druggability for {name} network genes...")
    druggability_results = {gene: assess_druggability(gene) for gene in centrality_df.index}
    centrality_df['is_druggable'] = [result[0] for result in druggability_results.values()]
    centrality_df['small_molecule'] = [result[1]['small_molecule'] for result in druggability_results.values()]
    centrality_df['antibody'] = [result[1]['antibody'] for result in druggability_results.values()]
    centrality_df['other_modality'] = [result[1]['other'] for result in druggability_results.values()]
    
    # Add gene information
    centrality_df['gene_symbol'] = centrality_df.index
    
    # Reorder columns
    column_order = ['gene_symbol', 'degree', 'betweenness', 'combined_score', 'is_hub', 
                    'is_druggable', 'small_molecule', 'antibody', 'other_modality']
    centrality_df = centrality_df[column_order]
    
    return centrality_df

# Analyze networks
up_centrality_df = analyze_network(G_up_targets, "up-regulated miRNA targets")
down_centrality_df = analyze_network(G_down_targets, "down-regulated miRNA targets")

# -----------------------------------------
# 6. Generate Interactive HTML Visualization
# -----------------------------------------
print("Generating interactive network visualizations...")

def visualize_network(G, centrality_df, name):
    if G.number_of_nodes() == 0:
        print(f"No nodes in {name} network. Skipping visualization.")
        return
    
    # Create a pyvis network
    net = Network(height="800px", width="100%", bgcolor="#222222", font_color="white")
    
    # Set physics layout
    net.barnes_hut(gravity=-80000, central_gravity=0.3, spring_length=250, spring_strength=0.001)
    
    # Add nodes to the network
    for node in G.nodes():
        # Get centrality and druggability information
        if node in centrality_df.index:
            info = centrality_df.loc[node]
            
            # Determine node size based on centrality (scaled for better visualization)
            size = 10 + (info['combined_score'] * 40)
            
            # Determine node color based on druggability
            if info['is_hub'] and info['is_druggable']:
                color = "#ff5722"  # Orange-red for druggable hubs
                title = f"HUB & DRUGGABLE: {node}<br>Centrality: {info['combined_score']:.3f}"
            elif info['is_hub']:
                color = "#e91e63"  # Pink for hubs
                title = f"HUB: {node}<br>Centrality: {info['combined_score']:.3f}"
            elif info['is_druggable']:
                color = "#4caf50"  # Green for druggable
                title = f"DRUGGABLE: {node}<br>Centrality: {info['combined_score']:.3f}"
            else:
                color = "#2196f3"  # Blue for regular nodes
                title = f"{node}<br>Centrality: {info['combined_score']:.3f}"
            
            # Add druggability details to the title
            if info['is_druggable']:
                modalities = []
                if info['small_molecule']: modalities.append("Small Molecule")
                if info['antibody']: modalities.append("Antibody")
                if info['other_modality']: modalities.append("Other Modality")
                title += f"<br>Druggable by: {', '.join(modalities)}"
        else:
            # Default values for nodes not in centrality_df
            size = 5
            color = "#9e9e9e"  # Grey
            title = node
        
        # Add the node
        net.add_node(node, title=title, size=size, color=color)
    
    # Add edges to the network
    for source, target, data in G.edges(data=True):
        # Edge width based on weight (if available)
        width = data.get('weight', 0.5) * 2
        net.add_edge(source, target, width=width, title=f"Confidence: {data.get('weight', 0.5):.2f}")
    
    # Set options
    net.set_options("""
    const options = {
      "nodes": {
        "borderWidth": 2,
        "borderWidthSelected": 4,
        "font": {
          "size": 12,
          "face": "Tahoma"
        }
      },
      "edges": {
        "color": {
          "inherit": true
        },
        "smooth": {
          "type": "continuous",
          "forceDirection": "none"
        }
      },
      "physics": {
        "forceAtlas2Based": {
          "gravitationalConstant": -50,
          "centralGravity": 0.01,
          "springLength": 100,
          "springConstant": 0.08
        },
        "minVelocity": 0.75,
        "solver": "forceAtlas2Based",
        "timestep": 0.5
      }
    }
    """)
    
    # Save the visualization
    output_path = os.path.join(network_dir, f"{name.replace(' ', '_')}_{timestamp}.html")
    net.save_graph(output_path)
    print(f"Interactive network visualization saved to: {output_path}")
    
    return output_path

# Generate visualizations
if len(up_centrality_df) > 0:
    up_viz_path = visualize_network(G_up_targets, up_centrality_df, "upregulated_mirna_targets")
else:
    up_viz_path = None
    
if len(down_centrality_df) > 0:
    down_viz_path = visualize_network(G_down_targets, down_centrality_df, "downregulated_mirna_targets")
else:
    down_viz_path = None

# -----------------------------------------
# 7. Export Results
# -----------------------------------------
print("Exporting results...")

# Export centrality and druggability data
if len(up_centrality_df) > 0:
    up_output_path = os.path.join(results_dir, f"upregulated_mirna_targets_analysis_{timestamp}.csv")
    up_centrality_df.to_csv(up_output_path, index=False)
    print(f"Upregulated miRNA targets analysis results saved to: {up_output_path}")

if len(down_centrality_df) > 0:
    down_output_path = os.path.join(results_dir, f"downregulated_mirna_targets_analysis_{timestamp}.csv")
    down_centrality_df.to_csv(down_output_path, index=False)
    print(f"Downregulated miRNA targets analysis results saved to: {down_output_path}")

# Export networks in GraphML format for use in Cytoscape/Gephi
if G_up_targets.number_of_nodes() > 0:
    graphml_path = os.path.join(network_dir, f"upregulated_mirna_targets_network_{timestamp}.graphml")
    nx.write_graphml(G_up_targets, graphml_path)
    print(f"Upregulated miRNA targets network saved in GraphML format to: {graphml_path}")

if G_down_targets.number_of_nodes() > 0:
    graphml_path = os.path.join(network_dir, f"downregulated_mirna_targets_network_{timestamp}.graphml")
    nx.write_graphml(G_down_targets, graphml_path)
    print(f"Downregulated miRNA targets network saved in GraphML format to: {graphml_path}")

# Create summary report
summary_file = os.path.join(results_dir, f"target_prioritization_summary_{timestamp}.txt")
with open(summary_file, 'w') as f:
    f.write("# miRNA Target Gene Prioritization Analysis\n\n")
    
    f.write("## Analysis Overview\n")
    f.write(f"- Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"- Total significant miRNAs: {len(df_mirna)}\n")
    f.write(f"- Upregulated miRNAs: {len(up_mirnas)}\n")
    f.write(f"- Downregulated miRNAs: {len(down_mirnas)}\n")
    f.write(f"- Target genes of upregulated miRNAs: {len(up_target_genes)}\n")
    f.write(f"- Target genes of downregulated miRNAs: {len(down_target_genes)}\n\n")
    
    f.write("## Network Analysis\n")
    f.write(f"- Upregulated miRNA targets network: {G_up_targets.number_of_nodes()} nodes, {G_up_targets.number_of_edges()} edges\n")
    f.write(f"- Downregulated miRNA targets network: {G_down_targets.number_of_nodes()} nodes, {G_down_targets.number_of_edges()} edges\n\n")
    
    f.write("## Hub Genes and Druggability\n")
    
    # Report on upregulated miRNA targets
    if len(up_centrality_df) > 0:
        hub_genes = up_centrality_df[up_centrality_df['is_hub']]
        druggable_hubs = hub_genes[hub_genes['is_druggable']]
        
        f.write(f"### Targets of Upregulated miRNAs (likely suppressed in HCC)\n")
        f.write(f"- Hub genes identified: {len(hub_genes)}\n")
        f.write(f"- Druggable hub genes: {len(druggable_hubs)}\n")
        
        # List top 10 druggable hub genes
        if len(druggable_hubs) > 0:
            f.write("\nTop druggable hub genes (by centrality score):\n")
            for i, (gene, row) in enumerate(druggable_hubs.head(10).iterrows(), 1):
                modalities = []
                if row['small_molecule']: modalities.append("Small Molecule")
                if row['antibody']: modalities.append("Antibody")
                if row['other_modality']: modalities.append("Other Modality")
                f.write(f"{i}. {row['gene_symbol']} (centrality: {row['combined_score']:.3f}, druggable by: {', '.join(modalities)})\n")
        f.write("\n")
    
    # Report on downregulated miRNA targets
    if len(down_centrality_df) > 0:
        hub_genes = down_centrality_df[down_centrality_df['is_hub']]
        druggable_hubs = hub_genes[hub_genes['is_druggable']]
        
        f.write(f"### Targets of Downregulated miRNAs (likely overexpressed in HCC)\n")
        f.write(f"- Hub genes identified: {len(hub_genes)}\n")
        f.write(f"- Druggable hub genes: {len(druggable_hubs)}\n")
        
        # List top 10 druggable hub genes
        if len(druggable_hubs) > 0:
            f.write("\nTop druggable hub genes (by centrality score):\n")
            for i, (gene, row) in enumerate(druggable_hubs.head(10).iterrows(), 1):
                modalities = []
                if row['small_molecule']: modalities.append("Small Molecule")
                if row['antibody']: modalities.append("Antibody")
                if row['other_modality']: modalities.append("Other Modality")
                f.write(f"{i}. {row['gene_symbol']} (centrality: {row['combined_score']:.3f}, druggable by: {', '.join(modalities)})\n")
        f.write("\n")
    
    f.write("## Interpretation\n")
    f.write("- Target genes of upregulated miRNAs are likely suppressed in HCC compared to normal tissue, suggesting potential tumor suppressor roles. Druggable hubs among these genes could be candidates for therapeutic enhancement or activation.\n")
    f.write("- Target genes of downregulated miRNAs are likely overexpressed in HCC compared to normal tissue, suggesting potential oncogenic roles. Druggable hubs among these genes are prime candidates for targeted inhibition.\n\n")
    
    f.write("## Output Files\n")
    if up_viz_path:
        f.write(f"- Interactive visualization of upregulated miRNA targets: {os.path.basename(up_viz_path)}\n")
    if down_viz_path:
        f.write(f"- Interactive visualization of downregulated miRNA targets: {os.path.basename(down_viz_path)}\n")
    
    if len(up_centrality_df) > 0:
        f.write(f"- Upregulated miRNA targets analysis results: {os.path.basename(up_output_path)}\n")
    if len(down_centrality_df) > 0:
        f.write(f"- Downregulated miRNA targets analysis results: {os.path.basename(down_output_path)}\n")
    
    f.write("\n## Advanced Visualization\n")
    f.write("For publication-quality visualizations with advanced edge bundling and layouts, the networks have been exported in GraphML format. These files can be imported into specialized network visualization tools:\n")
    f.write("- Cytoscape (https://cytoscape.org/): For interactive exploration and analysis\n")
    f.write("- Gephi (https://gephi.org/): For high-quality static visualizations\n")
    
    # Add explanatory note about the interactive visualizations
    f.write("\n## Note on Interactive Visualizations\n")
    f.write("The HTML visualizations provide a dynamic way to explore the networks:\n")
    f.write("- Node size represents centrality (larger = more central)\n")
    f.write("- Node colors represent: Orange-red (druggable hubs), Pink (non-druggable hubs), Green (druggable non-hubs), Blue (other genes)\n")
    f.write("- Edge thickness represents confidence of protein-protein interaction\n")
    f.write("- Hover over nodes and edges for detailed information\n")
    f.write("- Drag nodes to explore the network structure\n")
    f.write("- Use mouse wheel to zoom in/out\n")

print(f"Summary report saved to: {summary_file}")
print("Target prioritization analysis complete!") 