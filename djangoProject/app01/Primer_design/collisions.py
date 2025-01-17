from Bio.Seq import Seq
import Levenshtein

def hamming_distance(seq1, seq2):
    """Calculate the Hamming distance between two sequences."""
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def primer_payload_collision(primer, payloads, min_length=12, max_mismatches=2):
    """
    Check if the primer collides with any of the payloads.

    A collision happens if there are almost identical subsequences between the primer and payload
    that are longer than min_length and have at most max_mismatches.

    :param primer: The primer sequence (str).
    :param payloads: List of payload sequences (list of str).
    :param min_length: Minimum length of the subsequence to check (int).
    :param max_mismatches: Maximum number of mismatches allowed (int).
    :return: True if a collision is detected, False otherwise.
    """
    primer_len = len(primer)
    for payload in payloads:
        payload_len = len(payload)
        for i in range(primer_len - min_length + 1):
            for j in range(payload_len - min_length + 1):
                primer_subseq = primer[i:i + min_length]
                payload_subseq1 = payload[j:j + min_length - max_mismatches]
                if Levenshtein.distance(primer_subseq, payload_subseq1) <= max_mismatches:
                    if primer[i + min_length - max_mismatches:i + min_length] == payload[
                                                                                 j + min_length:j + min_length + max_mismatches]:
                        print(f"Collision detected with primer_subseq: {primer_subseq}")
                        print(f"Collision detected with payload: {payload[j:j + min_length]}")
                        return True
                payload_subseq2 = payload[j:j + min_length]
                if Levenshtein.distance(primer_subseq, payload_subseq2) <= max_mismatches:
                    print(f"Collision detected with primer_subseq: {primer_subseq}")
                    print(f"Collision detected with payload: {payload_subseq2}")
                    return True
    return False


def load_payloads_from_file(file_path):
    """
    Load payload sequences from a file.

    :param file_path: Path to the file containing payload sequences.
    :return: List of payload sequences (list of str).
    """
    with open(file_path, 'r') as f:
        payloads = f.read().splitlines()
    return payloads


def save_payloads_to_file(file_path, payloads):
    """
    Save modified payload sequences to a file.

    :param file_path: Path to the file to save modified payload sequences.
    :param payloads: List of modified payload sequences (list of str).
    """
    with open(file_path, 'w') as f:
        for payload in payloads:
            f.write(f"{payload}\n")


if __name__ == "__main__":
    forward_primer = "AAGGCAAGTTGTTACCAGCA"
    reverse_primer = "TGCGACCGTAATCAAACCAA"
    payload_file = './INNSE.dna'
    output_file = 'INNSE_Primer.dna'

    payloads = load_payloads_from_file(payload_file)

    # Convert the reverse primer to its reverse complement
    reverse_primer_rc = str(Seq(reverse_primer).reverse_complement())

    # Check for collisions
    print("Checking forward primer for collisions...")
    forward_collision = primer_payload_collision(forward_primer, payloads)

    print("Checking reverse primer for collisions...")
    reverse_collision = primer_payload_collision(reverse_primer_rc, payloads)

    if forward_collision or reverse_collision:
        print("Collision detected!")
    else:
        # Add primers to payloads
        modified_payloads = [forward_primer + payload + reverse_primer for payload in payloads]
        save_payloads_to_file(output_file, modified_payloads)
        print("Primers added to payloads and saved to output file.")
