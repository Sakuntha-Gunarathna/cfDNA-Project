# Author: Sakuntha Devaka Gunarathna
# Date: 2025-08-01


import pandas as pd
import os

# === User-defined settings ===
# Define start and end column names for the window (must match column headers in input files)
start_column_name = '<START_COLUMN>'  # e.g., '-500:-481'
end_column_name = '<END_COLUMN>'      # e.g., '480:499'

# Define output file name
output_file_name = 'Count_Matrix.txt'

# Get the current working directory (where input files are located)
dir_path = os.getcwd()
output_file_path = os.path.join(dir_path, output_file_name)

# Initialize an empty DataFrame for final results
final_df = pd.DataFrame()

# Process each file in the directory
for file in os.listdir(dir_path):
    # Skip unwanted files (scripts, text, bed, shell scripts, etc.). Please make sure only the Deeptools output files are in the folder.
    if file.endswith(('.py', '.txt', '.bed', '.sh')) or file in ['nohup.out', output_file_name]:
        continue

    print(f"Processing file: {file}")
    file_path = os.path.join(dir_path, file)

    try:
        # Read input data (skipping first 8 lines if they contain metadata)
        df = pd.read_csv(file_path, skiprows=8, header=0, delimiter='\t')
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        continue

    try:
        # Remove duplicate rows by Gene ID
        df = df.drop_duplicates(subset='Gene ID')

        # Select columns for the defined window
        columns_range = list(df.columns[df.columns.get_loc(start_column_name): df.columns.get_loc(end_column_name) + 1])
        
        # Calculate sum of signal across the window
        df['Row_Total'] = df[columns_range].sum(axis=1)

        # Derive sample name from file name (everything before first dot)
        sample_name = file.split('.')[0]

        # Initialize with Gene IDs if first file
        if final_df.empty:
            final_df['Gene ID'] = df['Gene ID']

        # Add counts for this sample
        final_df[sample_name] = df['Row_Total']

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        continue

# Save combined results
final_df.to_csv(output_file_path, sep='\t', index=False)
print(f"Data saved to {output_file_path}")
