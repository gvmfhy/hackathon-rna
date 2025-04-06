#!/usr/bin/env python
"""
Advanced Visualizations for miRNA Expression Analysis
Generates cutting-edge visualizations for miRNA expression data:
1. UMAP dimensionality reduction plot
2. miRNA-target network visualization
3. Interactive expression explorer
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx as nx
from sklearn.preprocessing import StandardScaler
import umap
import warnings
warnings.filterwarnings('ignore')

# Set output directory for results
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_dir = os.path.join(base_dir, "results")
os.makedirs(output_dir, exist_ok=True)

# Load the miRNA expression data
def load_data():
    print("Loading miRNA expression data...")
    tcga_file = os.path.join(base_dir, "data", "tcga", "TCGA_LIHC_miRNA_processed.csv")
    df = pd.read_csv(tcga_file)
    print(f"Dataset shape: {df.shape}")
    
    # Get counts of each group
    group_counts = df['Group'].value_counts()
    print("Group counts: ")
    print(group_counts)
    
    # Identify normal and tumor samples
    normal_samples = df[df['Group'] == 'Normal']
    tumor_samples = df[df['Group'] == 'NASH-HCC']
    
    print(f"Number of Normal samples: {len(normal_samples)}")
    print(f"Number of HCC samples: {len(tumor_samples)}")
    
    # Load DE results if available
    try:
        de_results_files = [f for f in os.listdir(output_dir) if f.startswith('mirna_de_significant') and f.endswith('.csv')]
        if de_results_files:
            # Get most recent file
            latest_file = sorted(de_results_files)[-1]
            significant_mirnas = pd.read_csv(os.path.join(output_dir, latest_file))
            print(f"Loaded {len(significant_mirnas)} significant miRNAs from {latest_file}")
        else:
            # Check for manually created file with names
            if os.path.exists(os.path.join(output_dir, "mirna_de_significant_with_names.csv")):
                significant_mirnas = pd.read_csv(os.path.join(output_dir, "mirna_de_significant_with_names.csv"))
                print(f"Loaded {len(significant_mirnas)} significant miRNAs from mirna_de_significant_with_names.csv")
            else:
                significant_mirnas = None
                print("No significant miRNA results file found.")
    except Exception as e:
        significant_mirnas = None
        print(f"Error loading DE results: {e}")
    
    # Load miRNA name mapping from the same place as mirna_differential_expression.py
    print("Loading miRNA name mapping...")
    mimat_to_name_map = {}
    
    try:
        # First try directly from mirbase dir
        mirbase_file = os.path.join(base_dir, "data", "mirbase", "hsa_miRBase.csv")
        if os.path.exists(mirbase_file):
            mirbase = pd.read_csv(mirbase_file)
            mimat_to_name_map = dict(zip(mirbase['MIMAT'], mirbase['Name']))
            print(f"Loaded {len(mimat_to_name_map)} mappings from miRBase file.")
        else:
            # If the file doesn't exist, try to read the mapping from a pre-processed file
            # that might have been created during differential expression analysis
            mimat_map_file = os.path.join(output_dir, "mimat_to_name_map.csv")
            if os.path.exists(mimat_map_file):
                mimat_map_df = pd.read_csv(mimat_map_file)
                mimat_to_name_map = dict(zip(mimat_map_df['MIMAT'], mimat_map_df['Name']))
                print(f"Loaded {len(mimat_to_name_map)} mappings from pre-processed mapping file.")
            else:
                # Try to find columns with names in the significant results
                if significant_mirnas is not None and 'Name' in significant_mirnas.columns and 'miRNA' in significant_mirnas.columns:
                    mimat_to_name_map = dict(zip(significant_mirnas['miRNA'], significant_mirnas['Name']))
                    print(f"Extracted {len(mimat_to_name_map)} mappings from significant results.")
                else:
                    print("No miRNA name mapping found. Will use MIMAT IDs as names.")
    except Exception as e:
        print(f"Error loading miRNA name mapping: {e}")
    
    # Load TargetScan data if available (but make it optional)
    targetscan = None
    try:
        targetscan_file = os.path.join(base_dir, "data", "targetscan", "TargetScanHuman.txt")
        
        if os.path.exists(targetscan_file):
            print("Loading TargetScan data...")
            # TargetScan has a lot of columns, we only need a few
            targetscan = pd.read_csv(targetscan_file, sep='\t', usecols=['miRNA', 'Gene Symbol', 'weighted context++ score'])
            print(f"Loaded TargetScan data with {len(targetscan)} miRNA-target pairs")
        elif os.path.exists(targetscan_file + ".zip"):
            # Unzip the file if needed
            print("Found zipped TargetScan data. Extracting...")
            import zipfile
            try:
                with zipfile.ZipFile(targetscan_file + ".zip", 'r') as zip_ref:
                    zip_ref.extractall(os.path.dirname(targetscan_file))
                    
                if os.path.exists(targetscan_file):
                    targetscan = pd.read_csv(targetscan_file, sep='\t', usecols=['miRNA', 'Gene Symbol', 'weighted context++ score'])
                    print(f"Loaded TargetScan data with {len(targetscan)} miRNA-target pairs")
            except zipfile.BadZipFile:
                print("Error: TargetScan zip file is corrupted or not a valid zip file.")
        else:
            print("TargetScan data not found. Network visualization will be skipped.")
    except Exception as e:
        print(f"Error loading TargetScan data: {e}")
        targetscan = None
    
    return df, normal_samples, tumor_samples, significant_mirnas, mimat_to_name_map, targetscan

def generate_umap_visualization(df, significant_mirnas=None, mimat_to_name_map=None):
    """Generate UMAP visualization to show sample clustering based on miRNA expression"""
    print("\nGenerating UMAP visualization...")
    
    # Extract expression data (all numeric columns)
    expression_data = df.select_dtypes(include=['float64', 'int64'])
    
    # If we have significant miRNAs, use only those for the visualization
    if significant_mirnas is not None and not significant_mirnas.empty:
        sig_mirnas = significant_mirnas['miRNA'].tolist()
        # Make sure all significant miRNAs are in the expression data
        sig_mirnas = [m for m in sig_mirnas if m in expression_data.columns]
        if sig_mirnas:
            print(f"Using {len(sig_mirnas)} significant miRNAs for UMAP.")
            expression_data = expression_data[sig_mirnas]
    
    # Scale the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(expression_data)
    
    # Apply UMAP
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
    umap_result = reducer.fit_transform(scaled_data)
    
    # Create a DataFrame with UMAP results
    umap_df = pd.DataFrame({
        'UMAP1': umap_result[:, 0],
        'UMAP2': umap_result[:, 1],
        'Group': df['Group'],
        'Sample_ID': df['Sample_ID'] if 'Sample_ID' in df.columns else df.index
    })
    
    # Create interactive UMAP plot with Plotly
    fig = px.scatter(
        umap_df, 
        x='UMAP1', 
        y='UMAP2', 
        color='Group',
        color_discrete_map={'Normal': 'blue', 'NASH-HCC': 'red'},
        hover_name='Sample_ID',
        title='UMAP Visualization of miRNA Expression Profiles',
        labels={'Group': 'Sample Type'},
        size_max=10,
        template='plotly_white'
    )
    
    # Add custom styling
    fig.update_traces(marker=dict(size=10, opacity=0.7, line=dict(width=1, color='white')))
    fig.update_layout(
        legend=dict(title='Sample Type', orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5),
        width=900,
        height=700,
        font=dict(family="Arial", size=14),
        hoverlabel=dict(font_size=12, font_family="Arial"),
        margin=dict(l=40, r=40, t=50, b=40)
    )
    
    # Save interactive UMAP plot
    umap_html_file = os.path.join(output_dir, f'umap_visualization_{timestamp}.html')
    fig.write_html(umap_html_file)
    print(f"Interactive UMAP visualization saved to: {umap_html_file}")
    
    # Also create a static version for publications
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='UMAP1', 
        y='UMAP2', 
        hue='Group', 
        data=umap_df,
        palette={'Normal': 'blue', 'NASH-HCC': 'red'},
        s=80,
        alpha=0.7
    )
    plt.title('miRNA Expression Profile Clustering', fontsize=16)
    plt.xlabel('UMAP Dimension 1', fontsize=14)
    plt.ylabel('UMAP Dimension 2', fontsize=14)
    plt.legend(title='Sample Type', fontsize=12, title_fontsize=14)
    
    # Customize the plot
    sns.despine()
    plt.tight_layout()
    
    # Save static plot
    umap_png_file = os.path.join(output_dir, f'umap_visualization_{timestamp}.png')
    plt.savefig(umap_png_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Static UMAP visualization saved to: {umap_png_file}")
    
    return umap_df

def generate_network_visualization(significant_mirnas, targetscan, mimat_to_name_map):
    """Generate network visualization showing miRNA-target interactions"""
    print("\nGenerating miRNA-target network visualization...")
    
    if significant_mirnas is None or targetscan is None:
        print("Cannot generate network visualization: missing miRNA or TargetScan data.")
        return
    
    # Get top miRNAs (most significant by p-value)
    top_mirnas = significant_mirnas.sort_values('padj').head(10)['miRNA'].tolist()
    
    # Convert MIMAT IDs to miRNA names if possible
    mirna_names = {}
    for mimat in top_mirnas:
        if mimat in mimat_to_name_map:
            # Clean the miRNA name (remove hsa- prefix for matching with TargetScan)
            mirna_name = mimat_to_name_map[mimat].replace('hsa-', '')
            mirna_names[mimat] = mirna_name
        else:
            mirna_names[mimat] = mimat
    
    # Build network of miRNA-target interactions
    G = nx.Graph()
    
    # Add miRNA nodes
    for mimat, name in mirna_names.items():
        G.add_node(name, type='miRNA', id=mimat)
    
    # Get targets for each miRNA and add edges
    for mimat, mirna_name in mirna_names.items():
        # Find targets for this miRNA
        # Note: TargetScan uses names without 'hsa-' prefix
        mirna_targets = targetscan[targetscan['miRNA'] == mirna_name]
        
        # Sort by interaction strength (lower score = stronger interaction)
        top_targets = mirna_targets.sort_values('weighted context++ score').head(5)
        
        # Add target nodes and edges
        for _, target in top_targets.iterrows():
            gene = target['Gene Symbol']
            score = target['weighted context++ score']
            
            # Add target node if not already in the graph
            if gene not in G:
                G.add_node(gene, type='target', id=gene)
            
            # Add edge with interaction strength
            # Convert score to edge weight (stronger interactions = thicker edges)
            # Scores are negative, where more negative = stronger interaction
            edge_weight = min(5, max(1, -score))
            G.add_edge(mirna_name, gene, weight=edge_weight, score=score)
    
    # If network is empty, return early
    if len(G.nodes()) == 0:
        print("No miRNA-target interactions found in the data.")
        return
    
    # Set node positions using a spring layout
    pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
    
    # Create the plot
    plt.figure(figsize=(14, 12))
    
    # Draw miRNA nodes (large red nodes)
    mirna_nodes = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'miRNA']
    nx.draw_networkx_nodes(G, pos, 
                          nodelist=mirna_nodes,
                          node_color='red',
                          node_size=800,
                          alpha=0.8)
    
    # Draw target nodes (smaller blue nodes)
    target_nodes = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'target']
    nx.draw_networkx_nodes(G, pos, 
                          nodelist=target_nodes,
                          node_color='skyblue',
                          node_size=400,
                          alpha=0.7)
    
    # Draw edges with varying thickness based on interaction strength
    for u, v, data in G.edges(data=True):
        width = data['weight']
        nx.draw_networkx_edges(G, pos, 
                              edgelist=[(u,v)],
                              width=width,
                              alpha=0.5,
                              edge_color='gray')
    
    # Add node labels
    nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold',
                           font_color='black')
    
    plt.title('miRNA-Target Interaction Network', fontsize=18)
    plt.text(0.01, 0.01, 'Edge thickness indicates interaction strength', 
             transform=plt.gca().transAxes, fontsize=12)
    
    # Create custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=15, label='miRNAs'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='skyblue', markersize=10, label='Target Genes')
    ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=12)
    
    # Remove axis
    plt.axis('off')
    plt.tight_layout()
    
    # Save the network visualization
    network_file = os.path.join(output_dir, f'mirna_target_network_{timestamp}.png')
    plt.savefig(network_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Network visualization saved to: {network_file}")

    # Create interactive version with Plotly
    try:
        edge_traces = []
        for edge in G.edges(data=True):
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            width = edge[2]['weight']
            
            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                line=dict(width=width, color='#888'),
                hoverinfo='none',
                mode='lines'
            )
            edge_traces.append(edge_trace)
        
        # Create node traces for miRNAs and targets
        mirna_trace = go.Scatter(
            x=[pos[node][0] for node in mirna_nodes],
            y=[pos[node][1] for node in mirna_nodes],
            text=[node for node in mirna_nodes],
            mode='markers+text',
            hoverinfo='text',
            marker=dict(
                color='red',
                size=20,
                line=dict(width=1, color='black')
            ),
            textposition="top center",
            name='miRNAs'
        )
        
        target_trace = go.Scatter(
            x=[pos[node][0] for node in target_nodes],
            y=[pos[node][1] for node in target_nodes],
            text=[node for node in target_nodes],
            mode='markers+text',
            hoverinfo='text',
            marker=dict(
                color='skyblue',
                size=15,
                line=dict(width=1, color='black')
            ),
            textposition="bottom center",
            name='Target Genes'
        )
        
        # Create the figure
        fig = go.Figure(data=edge_traces + [mirna_trace, target_trace],
                      layout=go.Layout(
                          title='miRNA-Target Interaction Network',
                          showlegend=True,
                          hovermode='closest',
                          margin=dict(b=0, l=0, r=0, t=40),
                          xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                          yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                          width=900,
                          height=900,
                          template='plotly_white'
                      ))
        
        # Add annotation explaining edge thickness
        fig.add_annotation(
            text="Edge thickness indicates interaction strength",
            xref="paper", yref="paper",
            x=0.01, y=0.01,
            showarrow=False
        )
        
        # Save interactive network
        network_html_file = os.path.join(output_dir, f'mirna_target_network_{timestamp}.html')
        fig.write_html(network_html_file)
        print(f"Interactive network visualization saved to: {network_html_file}")
    
    except Exception as e:
        print(f"Error creating interactive network visualization: {e}")

def generate_advanced_heatmap(df, significant_mirnas, mimat_to_name_map):
    """Generate an advanced interactive heatmap with clustering"""
    print("\nGenerating advanced interactive heatmap...")
    
    if significant_mirnas is None or significant_mirnas.empty:
        print("Cannot generate heatmap: no significant miRNAs found.")
        return
    
    # Select top miRNAs by p-value
    top_mirnas_ids = significant_mirnas.head(20)['miRNA'].tolist()
    
    # Get corresponding names, fallback to ID if name not found
    top_mirnas_names = [mimat_to_name_map.get(mid, mid) for mid in top_mirnas_ids]
    
    # LIMIT SAMPLES: Take fewer samples for a cleaner visualization
    num_normal_samples = min(15, len(df[df['Group'] == 'Normal']))
    num_tumor_samples = min(30, len(df[df['Group'] == 'NASH-HCC']))
    
    sampled_normal = df[df['Group'] == 'Normal'].sample(num_normal_samples, random_state=42) if len(df[df['Group'] == 'Normal']) > num_normal_samples else df[df['Group'] == 'Normal']
    sampled_tumor = df[df['Group'] == 'NASH-HCC'].sample(num_tumor_samples, random_state=42) if len(df[df['Group'] == 'NASH-HCC']) > num_tumor_samples else df[df['Group'] == 'NASH-HCC']
    plot_samples = pd.concat([sampled_normal, sampled_tumor])
    
    # Extract expression data
    heatmap_data = plot_samples[top_mirnas_ids].copy()
    
    # Z-score normalize the data
    from scipy.stats import zscore
    heatmap_z = pd.DataFrame(
        zscore(heatmap_data.values, axis=0),
        columns=top_mirnas_names,
        index=heatmap_data.index
    )
    
    # Create sample annotation
    sample_annotations = plot_samples['Group'].map({'Normal': 0, 'NASH-HCC': 1}).values
    
    # Create a clustered heatmap using hierarchical clustering
    from scipy.cluster.hierarchy import linkage, leaves_list
    
    # Cluster the samples (rows)
    row_linkage = linkage(heatmap_z.values, method='average', metric='euclidean')
    row_order = leaves_list(row_linkage)
    
    # Cluster the miRNAs (columns)
    col_linkage = linkage(heatmap_z.values.T, method='average', metric='euclidean')
    col_order = leaves_list(col_linkage)
    
    # Reorder the heatmap according to clustering
    heatmap_z_clustered = heatmap_z.iloc[row_order, col_order]
    
    # Reorder sample annotations too
    sample_annotations_clustered = [sample_annotations[i] for i in row_order]
    
    # Get ordered labels
    ordered_mirnas = [top_mirnas_names[i] for i in col_order]
    
    # Create interactive clustered heatmap with Plotly
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_z_clustered.values,
        x=ordered_mirnas,
        colorscale='RdBu_r',
        colorbar=dict(title='Z-score'),
        hovertemplate='miRNA: %{x}<br>Sample: %{y}<br>Z-score: %{z:.2f}<extra></extra>'
    ))
    
    # Add sample type annotation as a separate heatmap on the left
    fig.add_trace(go.Heatmap(
        z=[[a] for a in sample_annotations_clustered],
        colorscale=[[0, 'blue'], [1, 'red']],
        showscale=False,
        xaxis='x2'
    ))
    
    # Update layout to position the annotation heatmap
    fig.update_layout(
        title='Hierarchically Clustered miRNA Expression Patterns',
        xaxis=dict(
            domain=[0.05, 1],
            tickangle=45,
            tickfont=dict(size=12),
            title='miRNAs'
        ),
        xaxis2=dict(
            domain=[0, 0.05],
            showticklabels=False,
            title='Sample Type',
            anchor='y'
        ),
        yaxis=dict(
            showticklabels=False,
            title='Samples'
        ),
        width=1000,
        height=800,
        template='plotly_white',
        annotations=[
            dict(
                x=-0.02,
                y=0.5,
                xref='paper',
                yref='paper',
                text='Sample Type',
                showarrow=False,
                textangle=-90
            )
        ]
    )
    
    # Add a color legend for the sample types
    legend_annotations = [
        dict(
            x=-0.07,
            y=0.2,
            xref='paper',
            yref='paper',
            text='Normal',
            showarrow=False,
            font=dict(color='blue')
        ),
        dict(
            x=-0.07,
            y=0.1,
            xref='paper',
            yref='paper',
            text='HCC',
            showarrow=False,
            font=dict(color='red')
        )
    ]
    fig.update_layout(annotations=fig.layout.annotations + tuple(legend_annotations))
    
    # Save interactive heatmap
    heatmap_html_file = os.path.join(output_dir, f'advanced_heatmap_{timestamp}.html')
    fig.write_html(heatmap_html_file)
    print(f"Interactive heatmap saved to: {heatmap_html_file}")
    
def main():
    # Load all necessary data
    df, normal_samples, tumor_samples, significant_mirnas, mimat_to_name_map, targetscan = load_data()
    
    # If no significant miRNAs were found, create a minimal set for demonstration
    if significant_mirnas is None or len(significant_mirnas) == 0:
        print("\nNo significant miRNAs found. Creating a minimal set for visualization demo...")
        # Take top 10 miRNAs with highest variance across samples
        expression_data = df.select_dtypes(include=['float64', 'int64'])
        # Calculate variance for each miRNA
        var_series = expression_data.var()
        # Get top 10 miRNAs with highest variance
        top_vars = var_series.sort_values(ascending=False).head(10)
        
        # Create a minimal significant miRNAs dataframe
        significant_mirnas = pd.DataFrame({
            'miRNA': top_vars.index,
            'log2FC': [0.5] * len(top_vars),  # Dummy values
            'padj': [0.01] * len(top_vars)    # Dummy values
        })
        print(f"Created a demo set with {len(significant_mirnas)} high-variance miRNAs.")
    
    # Generate UMAP visualization
    try:
        umap_df = generate_umap_visualization(df, significant_mirnas, mimat_to_name_map)
        print("UMAP visualization completed successfully.")
    except Exception as e:
        print(f"Error generating UMAP visualization: {e}")
    
    # Generate miRNA-target network visualization only if TargetScan data is available
    if targetscan is not None and significant_mirnas is not None:
        try:
            generate_network_visualization(significant_mirnas, targetscan, mimat_to_name_map)
            print("Network visualization completed successfully.")
        except Exception as e:
            print(f"Error generating network visualization: {e}")
    else:
        print("Skipping network visualization due to missing data.")
    
    # Generate advanced interactive heatmap
    if significant_mirnas is not None:
        try:
            generate_advanced_heatmap(df, significant_mirnas, mimat_to_name_map)
            print("Interactive heatmap completed successfully.")
        except Exception as e:
            print(f"Error generating interactive heatmap: {e}")
    
    print("\nAdvanced visualizations complete!")

if __name__ == "__main__":
    main() 