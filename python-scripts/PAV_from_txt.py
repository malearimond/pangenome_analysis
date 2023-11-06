import os
import pandas as pd

# Define the directory containing your text files
text_files_dir = '/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/PAV/txt_files'

# Get a list of sample names from the text file names (excluding file extensions)
sample_names = [filename.split('.')[0] for filename in os.listdir(text_files_dir) if filename.endswith('.txt')]

# Initialize a dictionary to store gene names for each sample
sample_gene_sets = {sample: set() for sample in sample_names}

# Iterate through each text file and store gene names
for filename in os.listdir(text_files_dir):
    if filename.endswith('.txt'):
        sample_name = filename.split('.')[0]
        with open(os.path.join(text_files_dir, filename), 'r') as file:
            for line in file:
                sample_gene_sets[sample_name].add(line.strip())

# Create a set of all unique gene names across all samples
all_genes = set()
for gene_set in sample_gene_sets.values():
    all_genes.update(gene_set)

# Create a presence-absence matrix (list of lists)
presence_absence_matrix = []
for gene in all_genes:
    row = [1 if gene in sample_gene_sets[sample] else 0 for sample in sample_names]
    presence_absence_matrix.append(row)

# Create a DataFrame from the presence-absence matrix
df = pd.DataFrame(presence_absence_matrix, columns=sample_names)

# Save the DataFrame to a CSV file
output_file = 'presence_absence_matrix.csv'
df.to_csv(output_file, index=False)

print("Presence-absence matrix saved as {}".format(output_file))


