def remove_primers_from_payloads(payload_data, forward_primer_length, reverse_primer_length):
    """
    Remove forward and reverse primers from each sequence in the payload file.

    :param payload_file: Path to the file containing payload sequences.
    :param forward_primer_length: Length of the forward primer (int).
    :param reverse_primer_length: Length of the reverse primer (int).
    :param output_file: Path to the output file where modified sequences will be saved.
    """
    modified_sequences = []
    lines = payload_data.strip().split('\n')

    for line in lines:
        line = line.strip()
        # Remove forward primer
        payload_without_forward = line[forward_primer_length:]
        # Remove reverse primer
        payload_without_primers = payload_without_forward[:-reverse_primer_length]
        modified_sequences.append(payload_without_primers)

    return '\n'.join(modified_sequences)


if __name__ == "__main__":
    payload_file = './INNSE_Primer.dna'
    forward_primer_length = 20
    reverse_primer_length = 20

    with open(payload_file, 'r') as infile:
        payload_data = infile.read()

    modified_payloads = remove_primers_from_payloads(payload_data, forward_primer_length, reverse_primer_length)
