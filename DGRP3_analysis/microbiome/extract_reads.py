import os
import pandas as pd

# Define a function to extract data from a single file
def extract_data(file_path):
    data = {}
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
        sample_id = None
        for line in lines:
            if line.startswith("#SampleID"):
                sample_id = line.strip().split('\t')[1]
                data[sample_id] = {}
            elif line.startswith("k__") and "|s__" in line:
                parts = line.strip().split('\t')
                clade_name = parts[0]
                relative_abundance = parts[4] if len(parts) > 4 else 0.0
                data[sample_id][clade_name] = relative_abundance
    return data

# Define the directory where the files are located
directory = os.getcwd()

# Initialize an empty DataFrame to store the results
result_df = pd.DataFrame()

# Iterate through files in the directory
for subdir, _, files in os.walk(directory):
    for file in files:
        if file.endswith("_metaphlan_bugs_list.tsv"):
            file_path = os.path.join(subdir, file)
            print(file_path)
            sample_data = extract_data(file_path)
            sample_df = pd.DataFrame(sample_data).T
            result_df = pd.concat([result_df, sample_df], axis=0, sort=True)

# Fill missing values with 0
result_df.fillna(0, inplace=True)

# Transpose the DataFrame to have samples as columns
result_df = result_df.T

# Save the DataFrame to a CSV file
result_df.to_csv('clade_reads.csv', index=True)