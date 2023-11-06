def split_fasta(input_file, output_prefix, num_datasets):
    # Read the input FASTA file into a dictionary
    sequences = {}
    current_sequence = ''
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_sequence = line.strip()
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line.strip()

    # Shuffle the sequences randomly
    sequence_keys = list(sequences.keys())
    random.shuffle(sequence_keys)

    # Split the sequences into smaller datasets
    chunk_size = len(sequence_keys) // num_datasets
    for i in range(num_datasets):
        dataset_name = '{}_{}.fasta'.format(output_prefix, i)
        with open(dataset_name, 'w') as f:
            for key in sequence_keys[i * chunk_size:(i + 1) * chunk_size]:
                f.write(key + '\n')
                f.write(sequences[key] + '\n')

if __name__ == '__main__':
    input_file = 'input.fasta'  # Replace with your input FASTA file
    output_prefix = 'dataset'  # Prefix for output dataset filenames
    num_datasets = 5  # Number of smaller datasets to create

    split_fasta(input_file, output_prefix, num_datasets)
