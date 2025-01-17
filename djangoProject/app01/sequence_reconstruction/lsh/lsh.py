from app01.sequence_reconstruction.lsh.LSH_Clustering import *


def LSH_method(sequences):
    k_lsh = 3
    sim = 0.5  # number of MH signatures to be concatenated
    ell_lsh = int(1 / (sim ** k_lsh))  # number of LSH signatures
    m, k = 50, 7
    maxsig = 4 ** k

    print('****************************************Starting the main programme!***************************************')
    # original_sequences, reads = read_file(encoding_file, sequencing_file, min_length=100, max_length=120)

    data_list = processing_data(sequences)

    # s1 = time.time()
    minhash = Minhash_sign(m, k)
    kmer_indexes = []
    minhash_signs = []
    for seq in data_list:
        kmer_index, minhash_sign = minhash.generate_signature(seq[:40])
        kmer_indexes.append(kmer_index)
        minhash_signs.append(minhash_sign)
    pairs = extract_similar_pairs(minhash_signs, m, k_lsh, ell_lsh, maxsig)
    clusters = center_cluster(pairs)
    final_clusters = [[c] + list(set(clusters[c])) for c in clusters if len(clusters[c]) >= 3]
    th = 35  # filtering threshold
    fclusts = []
    for i, c in enumerate(final_clusters):
        cent = data_list[c[0]]
        cc = [c[0]]
        for e in c[1:]:
            score = max_match(cent, data_list[e])
            if score >= th:
                cc += [e]
        fclusts += [cc]
        if i % 1000 == 0:
            print("%", round(i * 100 / len(final_clusters), 2), "of the clusters are filtered.")

    output_list_noseq = ""
    for i, cluster in enumerate(fclusts):
        output_list_noseq += f">Cluster{i + 1}\r\n"
        for index in cluster:
            output_list_noseq += f"{data_list[index]}\r\n"
    output_list = ""
    sequences_index_dict = {}
    i = 0
    for line in sequences:
        line = line.rstrip()
        if line.startswith('>'):
            seq_label = int(re.search(r'Seq(\d+)', line).group(1))
            # print(i, seq_label)
            sequences_index_dict[i] = seq_label
            i += 1
        else:
            continue
    for i, cluster in enumerate(final_clusters):
        for seq in cluster:
            output_list += f'>Seq{sequences_index_dict[seq]}_Cluster{i + 1}\r\n'
            output_list += f'{data_list[seq]}\r\n'


    return output_list_noseq, output_list

