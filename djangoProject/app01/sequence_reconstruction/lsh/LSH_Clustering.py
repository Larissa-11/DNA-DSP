import numpy as np
from random import shuffle
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA
import itertools
import operator
import time
import zipfile
import subprocess as sp
import hashlib
# from include import *
from random import random
import operator
import math
from sklearn import metrics
import re


# =====Convert fasta files to lists=====
def fasta_to_list(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.rstrip()
            sequence = ''
            if line.startswith('>'):
                continue
            else:
                sequence += line
            sequences.append(sequence)
    return sequences


# =====Read and process sequencing files=====
def processing_data(line_list, min_length=100, max_length=300):
    # original_sequences = fasta_to_list(original_file)
    sequencing_sequences = []
    for line in line_list:
        line = line.rstrip()
        sequence = ''
        if line.startswith('>'):
            continue
        else:
            sequence += line
        if min_length <= len(sequence) <= max_length:
            sequencing_sequences.append(sequence)
    return sequencing_sequences


# =====Perform k-mer slicing of the sequence and compute the k-mer index of the sequence=====
def DNA_kmer_index(seq, k=3):  # The default step size is 1
    kmer = []
    for ell in range(len(seq) - k + 1):
        nstr = seq[ell:ell + k]
        index = 0
        for j, c in enumerate(nstr):
            if c == 'A':
                i = 0
            elif c == 'C':
                i = 1
            elif c == 'G':
                i = 2
            elif c == 'T':
                i = 3
            else:
                index = -1
                break
            index += i * (4 ** j)
        kmer += [index]
    return kmer


# =====min-hash object=====#
class Minhash_sign():
    # min-hash of k-mers
    def __init__(self, m, k):
        # table is a random hash table
        # m is the number of signatures
        self.tables = [np.random.permutation(4 ** k) for i in range(m)]
        self.k = k

    # Generate a list of minimum hash signatures
    def generate_signature(self, seq):
        kmer_index = DNA_kmer_index(seq, self.k)
        minhash_sign = [min([table[i] for i in kmer_index]) for table in self.tables]
        return kmer_index, minhash_sign


# =====Calculate the LSH signatures of the sequences based on the list of minimum hash signatures of the sequences =====
# =====and extract pairs of sequences with equal LSH signatures=====
def extract_similar_pairs(sigs, m, k_lsh, ell_lsh, maxsig):
    # sigs: minhash signatures 最小哈希签名列表
    # m: 最小哈希签名列表中的签名数
    # ell_lsh: number of LSH signatures 要生成的 LSH 签名数
    # k_lsh: number of MH signatures to be concatenated 每个 LSH 签名中要连接的最小哈希签名数
    # we use generatrs to yield a number of pairs at a time for the sake of memory efficiency
    pairs = set([])
    # generate ell_lsh random indices
    for ell in range(ell_lsh):
        pair_count = 0
        s = time.time()
        lshinds = np.random.permutation(m)[:k_lsh]
        # generate LSh signatures
        lshsigs = []
        for sig in sigs:
            lshsig = 0
            for i, lshind in enumerate(lshinds):
                lshsig += sig[lshind] * (maxsig ** i)
            lshsigs += [lshsig]
        d = {}
        for ind, sig in enumerate(lshsigs):
            if sig in d:
                d[sig] += [ind]
            else:
                d[sig] = [ind]
        for candidates in d.values():
            cent = set([])
            if len(candidates) > 1:
                for pair in itertools.combinations(candidates, 2):
                    cent.add(pair[0])
                    if len(cent) == 1:
                        pairs.add(pair)
                    else:
                        break
        yield pairs, ell
        pair_count += len(pairs)
        pairs = set([])


# =====Generate clusters based on extracted sequence pairs=====#
def center_cluster(pairs):
    clusters = {}
    pairsize = 0
    t_counter = 0
    hold = 0
    ell_copy = 0
    while not hold:
        try:
            out = next(pairs)
            # print(out)
            pairs_sort = list(out[0])
            ell = out[1]
            pairsize += len(pairs_sort)
            pairs_sort.sort()
            star_time = time.time()
            for (u, v) in pairs_sort:
                if u in clusters:
                    clusters[u] += [v]
                if v in clusters:
                    clusters[v] += [u]
                if v not in clusters and u not in clusters:
                    clusters[u] = [v]
                #
                # if u not in clusters:
                #     clusters[u] = [v]
                # else:
                #     clusters[u] += [v]
            t_counter += time.time() - star_time
            ell_copy = ell
            print("Clustering time for LSH", ell_copy, ":", t_counter)
        except StopIteration:
            hold = 1
            print("clustering completed", "---", pairsize, "pairs clustered")
    return clusters


# =====max matching=====#
def max_match(seq1, seq2):
    # This function checks whether seq1 and seq2 are similar or not
    # Checking all pairs within a cluster dramatically increases the time complexity,
    # so by default, in the next cell, we call this function to only check the pairs
    # that one of their members is the cluster center

    alignment, score, start_end_positions \
        = local_pairwise_align_ssw(DNA(seq1), DNA(seq2), match_score=2, mismatch_score=-3)
    a = str(alignment[0])
    b = str(alignment[1])
    ctr = 0
    for i, j in zip(a, b):
        if i == j:
            ctr += 1
    return ctr



