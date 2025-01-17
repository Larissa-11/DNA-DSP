from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import random


def gc_content(seq):
    """Calculate the GC content of the sequence."""
    return (seq.count('G') + seq.count('C')) / len(seq)


def has_homopolymers(seq, max_len):
    """Check if the sequence contains homopolymers longer than max_len."""
    for base in 'ACGT':
        if base * (max_len + 1) in seq:
            return True
    return False


def hamming_distance(seq1, seq2):
    """Calculate the Hamming distance between two sequences."""
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def reverse_complement(seq):
    """Generate the reverse complement of a sequence."""
    return str(Seq(seq).reverse_complement())


def check_inter_complementarity(new_seq, primer_library, inter_complementarity):
    """Check for inter-complementarity between a new sequence and existing primers."""
    new_rev_comp = reverse_complement(new_seq)
    for existing_seq in primer_library:
        for i in range(len(new_seq) - inter_complementarity + 1):
            for j in range(len(existing_seq) - inter_complementarity + 1):
                if new_seq[i:i + inter_complementarity] == reverse_complement(existing_seq[j:j + inter_complementarity]):
                    return True
    return False


def check_intra_complementarity(seq, intra_complementarity):
    """Check for intra-complementarity within a sequence."""
    rev_comp = reverse_complement(seq)
    for i in range(len(seq) - intra_complementarity + 1):
        if seq[i:i + intra_complementarity] in rev_comp:
            return True
    return False


def generate_sequence(length):
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def design_primers(params):
    """Design primers that satisfy all the given constraints."""
    primers = []
    # i = 1
    while True:
        # print("运行次数：", i)
        # i+=1
        seq = generate_sequence(params['length'])

        # Check GC content
        a = params['gc_content_min']
        b = gc_content(seq)
        c = params['gc_content_max']
        if not (params['gc_content_min'] <= gc_content(seq) <= params['gc_content_max']):
            continue

        # Check melting temperature
        temp = mt.Tm_NN(seq)
        if not (params['Tm_min'] <= temp <= params['Tm_max']):
            continue

        # Check homopolymers
        if has_homopolymers(seq, params['homo_max_len']):
            continue

        # Check Hamming distance with existing primers
        if any(hamming_distance(seq, primer) < params['hamming_distance'] for primer in primers):
            continue

        # Check inter- and intra-complementarity
        if check_intra_complementarity(seq, params['intra_complementarity']) or check_inter_complementarity(seq, primers, params['inter_complementarity']):
            continue

        primers.append(seq)
        # Stop when we have at least one primer
        if len(primers) >= params['number']:
            break
    return primers

# Main function to get user input and design primers
if __name__ == "__main__":
    params = {
        'length': 20,
        'gc_content_min': 0.45,
        'gc_content_max': 0.55,
        'homo_max_len': 3,
        'Tm_min': 55,
        'Tm_max': 60,
        'hamming_distance': 6,
        'inter_complementarity': 10,
        'intra_complementarity': 4,
        'number': 5
    }

    primers = design_primers(params)
    for i, primer in enumerate(primers):
        print(f"Primer {i + 1}: {primer}")
