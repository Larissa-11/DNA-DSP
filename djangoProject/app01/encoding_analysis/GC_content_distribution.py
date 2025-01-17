import random
import pandas as pd

from app01.encoding_analysis import data_handle

class GcContentDistribution:
    def compute_gc_content(self, sequence):
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = gc_count / len(sequence) * 100
        return gc_content

    def gc_content_distribution(self, dna_sequences_list):
        gc_contents_all = []
        for dna_sequences in dna_sequences_list:
            gc_contents = [self.compute_gc_content(sequence) for sequence in dna_sequences]
            gc_contents_all.append(gc_contents)
        return gc_contents_all

    def local_gc_content_distribution(self, dna_sequences_list):
        # 选择随机 DNA 片段并计算不同长度的 GC 含量
        sequence_lengths = list(range(6, 151, 6))

        gc_data = {'Sequence_Length': [], 'GC_Content': []}
        gc_data['Sequence_Length'].extend(sequence_lengths)
        for dna_sequences in dna_sequences_list:
            random_sequence = random.choice(dna_sequences)

            for length in sequence_lengths:
                sequence_subset = random_sequence[:length]
                gc_content = self.compute_gc_content(sequence_subset)
                gc_data['GC_Content'].append(gc_content)
        return gc_data
