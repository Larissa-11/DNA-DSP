from numpy import fromfile, uint8

from app01.encoding_analysis import data_handle


def calculate_net_information_density(bit_size, dna_sequences_list):
    net_info_densities = []
    for dna_sequences in dna_sequences_list:
        a = len(dna_sequences)
        b = len(dna_sequences[0])
        nt_size = a * b
        net_info_density = bit_size / nt_size
        net_info_density_rounded = round(net_info_density, 2)
        net_info_densities.append(net_info_density_rounded)

    return net_info_densities
def find_longest_homopolymer(dna_sequences_list):
    max_lengths = []
    for sequences in dna_sequences_list:
        max_length = 1
        for sequence in sequences:
            current_length = 1
            for i in range(1, len(sequence)):
                if sequence[i] == sequence[i - 1]:
                    current_length += 1
                    if current_length > max_length:
                        max_length = current_length
                else:
                    current_length = 1
        max_lengths.append(max_length)
    return max_lengths
def motifs(dna_sequences_list):
    motifs_counts = []
    for sequences in dna_sequences_list:
        motifs_count = 0
        for sequence in sequences:
            motifs = ["GGC", "GAATTC"]
            for missing_segment in motifs:
                if missing_segment in "".join(sequence):
                    motifs_count = "".join(sequence).count(missing_segment) + motifs_count
        motifs_counts.append(motifs_count)
    return motifs_counts
def read_bits_from_file(path, need_logs=False):
    if need_logs:
        print("Read binary matrix from file: " + path)
    values = fromfile(file=path, dtype=uint8)
    if need_logs:
        print("There are " + str(len(values) * 8) + " bits in the inputted file. ")
    return len(values) * 8


if __name__ == "__main__":
    img_path = './0801.png'
    DNA_paths = ['./0801.dna', './DNA-QLC.dna', './DNAFountain.dna']
    with open(img_path, mode='rb') as f:
         bit_size = read_bits_from_file(f)
    dna_sequences = [data_handle.read_dna_file(DNA_path) for DNA_path in DNA_paths]
    net_info_density = calculate_net_information_density(bit_size, dna_sequences)
    print("Net information density for all files: ", net_info_density)
    max_lengths = find_longest_homopolymer(DNA_paths)
    print("The length of the homopolymer：", max_lengths)
    motifs_count = motifs(DNA_paths)
    print('The motifs_count is :', motifs_count)
