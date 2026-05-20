import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances
import umap
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

def load_data(input_file):
    """
    Load input CSV file containing microbiome abundance data.
    Assumes first 7 columns contain taxonomic names and columns 8 onwards are relative abundances.
    """
    data = pd.read_csv(input_file)
    taxa = data.iloc[:, :7]
    abundances = data.iloc[:, 7:]
    return taxa, abundances

def perform_pca(data):
    """
    Perform PCA on the given data.
    """
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data.T)  # Transpose for samples in rows
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)
    explained_variance = pca.explained_variance_ratio_
    return pca_result, explained_variance

def perform_mds(data):
    """
    Perform MDS using Bray-Curtis distance, handling NaNs in the distance matrix.
    """
    # Transpose data to have samples in rows
    data = data.T

    # Remove samples with zero sum
    valid_samples = data.sum(axis=1) > 0
    data = data[valid_samples]
    if data.shape[0] == 0:
        raise ValueError("All samples have zero total abundance. Cannot compute Bray-Curtis distances.")

    # Compute Bray-Curtis distance
    bray_curtis_dist = pairwise_distances(data, metric="braycurtis")

    # Replace NaNs with maximum distance (1.0 for Bray-Curtis)
    bray_curtis_dist = np.nan_to_num(bray_curtis_dist, nan=1.0)

    # Perform MDS
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
    mds_result = mds.fit_transform(bray_curtis_dist)
    return mds_result

def perform_umap(data):
    """
    Perform UMAP on the given data.
    """
    umap_model = umap.UMAP(n_components=2, random_state=42)
    umap_result = umap_model.fit_transform(data.T)  # Transpose for samples in rows
    return umap_result

def plot_embedding(embedding, method, output_file, explained_variance=None):
    """
    Plot a 2D embedding (PCA, MDS, or UMAP) and save it to a file.
    """
    embedding_df = pd.DataFrame(embedding, columns=["Component 1", "Component 2"])
    embedding_df["Sample"] = embedding_df.index.astype(str)

    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=embedding_df,
        x="Component 1",
        y="Component 2",
        s=50,
        alpha=0.7
    )
    plt.title(f"{method} Embedding", fontweight="bold")
    plt.xlabel("Component 1")
    plt.ylabel("Component 2")
    if explained_variance is not None:
        plt.xlabel(f"Component 1 ({explained_variance[0]*100:.2f}%)")
        plt.ylabel(f"Component 2 ({explained_variance[1]*100:.2f}%)")
    plt.savefig(output_file, dpi=300)
    plt.close()

def main(input_file, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    taxa, abundances = load_data(input_file)

    # Perform PCA
    pca_result, explained_variance = perform_pca(abundances)
    plot_embedding(pca_result, "PCA", os.path.join(output_dir, "pca_plot.png"), explained_variance)

    # Perform MDS
    mds_result = perform_mds(abundances)
    plot_embedding(mds_result, "MDS (Bray-Curtis)", os.path.join(output_dir, "mds_plot.png"))

    # Perform UMAP
    umap_result = perform_umap(abundances)
    plot_embedding(umap_result, "UMAP", os.path.join(output_dir, "umap_plot.png"))

    print(f"Results saved in {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform PCA, MDS, and UMAP on microbiome relative abundance data.")
    parser.add_argument("-I", "--input", type=str, required=True, help="Input CSV file path")
    parser.add_argument("-O", "--output", type=str, required=True, help="Output directory for plots")
    
    args = parser.parse_args()

    # Check if input file exists
    if not os.path.isfile(args.input):
        print(f"Error: The input file '{args.input}' does not exist.")
        exit(1)

    main(args.input, args.output)
