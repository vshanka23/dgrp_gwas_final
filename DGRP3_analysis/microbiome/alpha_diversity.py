import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import entropy
from skbio.diversity.alpha import chao1
from scipy.special import gammaln
import argparse
import os

# Function to calculate Shannon diversity
def shannon_diversity(counts):
    return entropy(counts, base=2)

# Function to calculate Simpson diversity
def simpson_diversity(counts):
    p = counts / counts.sum()
    return 1 - np.sum(p ** 2)

# Function to calculate Simpson's Evenness (E)
def simpsons_evenness(counts):
    simpson = simpson_diversity(counts)
    richness = np.count_nonzero(counts)
    return simpson / (richness - 1) if richness > 1 else 0

# Function to calculate McIntosh's diversity index
def mcintosh_diversity(counts):
    counts = np.maximum(counts, 0)
    S = np.count_nonzero(counts)
    p = counts / counts.sum()
    return (S - 1) / (S - 1 + np.sum(p ** 2)) if S > 1 else 0

# Function to calculate Gini index
def gini_index(counts):
    counts = np.sort(counts)
    n = len(counts)
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * counts) - (n + 1) * np.sum(counts)) / (n * np.sum(counts))

# Function to calculate Pielou's evenness
def pielou_evenness(counts):
    shannon = shannon_diversity(counts)
    richness = np.count_nonzero(counts)
    return shannon / np.log(richness) if richness > 0 else 0

def renyi_entropy(counts, alpha):
    counts = np.maximum(counts, 0)
    p = counts / counts.sum()
    if alpha == 1:
        return shannon_diversity(counts)
    return (1 / (1 - alpha)) * np.log(np.sum(p ** alpha))

# Function to calculate Brillouin's diversity index
def brillouin_diversity(counts):
    counts = np.maximum(counts, 0)
    N = np.sum(counts)
    if N == 0:
        return 0
    S = len(counts)
    H_B = gammaln(N + 1) - np.sum(gammaln(counts + 1)) - np.log(N)
    return H_B

# Function to calculate Chao1 index
def chao1_index(counts):
    return chao1(counts)

def main(input_file, output_plot, output_csv):
    # Load the data
    data = pd.read_csv(input_file)

    # Extract the relative abundances (columns 8 onwards)
    abundances = data.iloc[:, 7:]

    # Calculate alpha diversity metrics for each sample
    alpha_diversity_metrics = pd.DataFrame({
        "Sample": abundances.columns,
        "Shannon's D": abundances.apply(shannon_diversity, axis=0),
        "Simpson Index": abundances.apply(simpson_diversity, axis=0),
        "Simpson's E": abundances.apply(simpsons_evenness, axis=0),
        "McIntosh's D": abundances.apply(mcintosh_diversity, axis=0),
        "Gini Index": abundances.apply(gini_index, axis=0),
        "Pielou's Evenness": abundances.apply(pielou_evenness, axis=0),
        "Chao1 Index": abundances.apply(chao1_index, axis=0),
        "Renyi's Entropy ": abundances.apply(lambda x: renyi_entropy(x, 2), axis=0),
        "Brillouin Index": abundances.apply(brillouin_diversity, axis=0)
    })

    # Save numerical values to a CSV file
    alpha_diversity_metrics.to_csv(output_csv, index=False)

    # Reshape the data for plotting
    melted_metrics = alpha_diversity_metrics.melt(id_vars='Sample', var_name='Metric', value_name='Value')

    # Generate multipanel violin plots with separate scales for each panel
    g = sns.FacetGrid(melted_metrics, col="Metric", col_wrap=3, sharex=False, sharey=False)
    g.map(sns.violinplot, "Value", inner="quartile", scale="width")
    g.set_titles(col_template="{col_name}", fontweight='bold')
    
    for ax in g.axes.flat:
        ax.set_xlabel("Index Value")
        ax.set_ylabel("Density")

    # Add a title
    plt.subplots_adjust(left=0.1, top=0.9, hspace=0.5, wspace=0.4)
    g.fig.suptitle('Measures of microbiome diversity and richness', fontweight='bold')

    # Save the plot to the output file with high resolution
    plt.savefig(output_plot, dpi=600)
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate alpha diversity metrics and generate multipanel violin plots.')
    parser.add_argument('-I', '--input', type=str, required=True, help='Input CSV file path')
    parser.add_argument('-O', '--output', type=str, required=True, help='Output file path for the plot')
    parser.add_argument('-C', '--csv_output', type=str, required=True, help='Output CSV file path for numerical values')

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.isfile(args.input):
        print(f"Error: The input file '{args.input}' does not exist.")
        exit(1)

    main(args.input, args.output, args.csv_output)