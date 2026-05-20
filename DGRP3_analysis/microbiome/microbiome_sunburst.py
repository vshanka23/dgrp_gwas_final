import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import sys

def create_core_microbiome_sunburst(abundance_table_path, output_plot_path):
    # Load the abundance table
    feature_table = pd.read_csv(abundance_table_path)
    
    # Removing viruses to make automatic color plotting easier
    #feature_table = feature_table[feature_table['Kingdom']!="Viruses"]

    # Identify taxonomic levels and abundance columns dynamically
    taxa_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    abundance_columns = feature_table.columns[len(taxa_levels):]

    # Prefix taxonomic level names with appropriate initials
    for level in taxa_levels:
        feature_table[level] = level[:1].lower() + '_' + feature_table[level].astype(str)

    # Calculate presence/absence and filter core taxa
    threshold = 0.00
    presence_absence = feature_table[abundance_columns] > 0
    core_taxa_mask = presence_absence.mean(axis=1) >= threshold
    core_feature_table = feature_table[core_taxa_mask]

    # Compute average abundance across all samples
    core_feature_table['AverageAbundance'] = core_feature_table[abundance_columns].mean(axis=1)

    # Summarize data at each taxonomic level for the core microbiome
    core_summary = core_feature_table.groupby(taxa_levels)['AverageAbundance'].sum().reset_index()

    # Create a dictionary to hold the hierarchical data
    hierarchical_data = {}

    def add_to_hierarchy(hierarchy, levels, abundance):
        current_level = hierarchy
        for level in levels[:-1]:  # Traverse until the second last level
            if level not in current_level:
                current_level[level] = {'children': {}, 'abundance': 0}
            current_level[level]['abundance'] += abundance
            current_level = current_level[level]['children']
        # The last level should not have 'children', only 'abundance'
        final_level = levels.iloc[-1]
        if final_level not in current_level:
            current_level[final_level] = {'abundance': 0}
        current_level[final_level]['abundance'] += abundance

    for _, row in core_summary.iterrows():
        add_to_hierarchy(hierarchical_data, row[taxa_levels], row['AverageAbundance'])

    # Convert hierarchical data to a format suitable for plotting
    def hierarchy_to_sunburst(hierarchy, parent_name=''):
        labels = []
        parents = []
        values = []
        def recurse(hierarchy, parent_name):
            for name, data in hierarchy.items():
                labels.append(name)
                parents.append(parent_name)
                values.append(data['abundance'])
                if 'children' in data and data['children']:
                    recurse(data['children'], name)
        recurse(hierarchy, parent_name)
        return labels, parents, values

    labels, parents, values = hierarchy_to_sunburst(hierarchical_data)

    # Create sunburst plot
    fig = go.Figure(go.Sunburst(
        labels=labels,
        parents=parents,
        values=values,
        branchvalues='total',
        maxdepth=6  # Adjust this to show desired depth
    ))

    # Update layout
    fig.update_layout(margin=dict(t=0, l=0, r=0, b=0)) 
    #,sunburstcolorway = px.colors.qualitative.Pastel

    #new_labels = [label.split('_', 1)[-1] for label in labels]
    #fig.update_traces(labels=new_labels)
    

    # Export the figure as PNG with 300 DPI
    fig.write_image(output_plot_path, format='svg', scale=6)  # scale parameter to achieve 300 DPI
    
    print(f"Sunburst plot saved as {output_plot_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python microbiome_sunburst.py <path to abundance table> <path and name of the sunburst plot>")
        sys.exit(1)
    
    abundance_table_path = sys.argv[1]
    output_plot_path = sys.argv[2]
    create_core_microbiome_sunburst(abundance_table_path, output_plot_path)