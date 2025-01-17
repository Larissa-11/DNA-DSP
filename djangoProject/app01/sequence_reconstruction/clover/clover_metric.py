import time
from src.compute_performances import *
from src.sequencing_preprocessing import *
import operator
import subprocess
from random import shuffle
from Bio import AlignIO
from tqdm import tqdm


def read_cluster_pair(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    # 使用eval函数将字符串转换为元组列表
    original_data = eval(content)
    # 将字符串中的数字字符串转换为整数
    converted_data = [(int(x), int(y)) for x, y in original_data]
    return converted_data


def center_cluster(pairs):
    clusters = {}
    for (u, v) in pairs:
        if u in clusters:
            clusters[u] += [v]

        if v in clusters:
            clusters[v] += [u]

        if v not in clusters and u not in clusters:
            clusters[v] = [u]

    return clusters


def muscle(cluster):
    input_alignment = "output/msa/clm.fasta"
    output_alignment = "output/msa/clmout.fasta"
    # write cluster to file
    file = open(input_alignment, "w")
    for i, c in enumerate(cluster):
        file.write(">S%d\n" % i)
        file.write(c)
        file.write("\n")
    file.close()
    # build and execute the muscle5 command
    muscle_exe = "muscle5"
    muscle_command = [muscle_exe, "-align", input_alignment, "-output", output_alignment]
    # muscle_command = [muscle_exe, "-super5", input_alignment, "-output", output_alignment]
    try:
        subprocess.run(muscle_command, check=True)
        print("MUSCLE alignment completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE: {e}")
    # read the aligned sequences
    msa_sequences = AlignIO.read(output_alignment, "fasta")
    aligned_cluster = []
    for i in msa_sequences:
        aligned_cluster += [i.seq]
    return aligned_cluster


def clusters_msa(clusters, reads, maxsize=15):
    # align clusters, generate candidates
    msa_results = []
    for i, cluster_indexes in enumerate(tqdm(clusters, ncols=100)):
        cluster = [reads[index] for index in cluster_indexes]
        if len(cluster) < 3:
            continue
        elif len(cluster) <= maxsize:
            aligned_cluster = muscle(cluster)
            msa_results.append(aligned_cluster)
        else:
            for j in range(5):
                shuffle(cluster)
                aligned_cluster = muscle(cluster[:maxsize])
                msa_results.append(aligned_cluster)
        # Alignment of progress
        if i % 1000 == 0:
            print("%", round(i * 100 / len(clusters), 2), "of the clusters are aligned.")
    return msa_results


def sequence_voting(msa, weight=0.4):
    # assume reads have the same length
    res = ""
    for i in range(len(msa[0])):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0, 'N': 0}
        for j in range(len(msa)):
            counts[msa[j][i]] += 1
        counts['-'] *= weight
        mv = max(counts.items(), key=operator.itemgetter(1))[0]
        if mv != '-':
            res += mv
    return res


def generate_candidates(msa_results):
    candidates = []
    for msa in msa_results:
        candidates.append(sequence_voting(msa, weight=0.5))
    return candidates


if __name__ == '__main__':
    iteration = 1
    # total_fraction = []
    # total_re_rate = []
    for it in range(iteration):
        # encoding_file = '../datasets/test_data/test_origin_50.fasta'
        # sequencing_file = '../datasets/test_data/test_cluster_50_with_labels.fasta'
        # encoding_file = '../datasets/encoding/nanopore_centers.fasta'
        # sequencing_file = '../datasets/clustering/nanopore_clusters_with_labels.fasta'
        encoding_file = '../datasets/encoding/P10_5_BDDP210000009_1A_encoding.fasta'
        sequencing_file = '../datasets/sequencing/P10_5_BDDP210000009_1A_sequencing_with_labels.fasta'

        correct_length = 200
        cluster_file = 'output/clustering_output.txt'
        original_sequences, reads = read_file(encoding_file, sequencing_file, correct_length)
        # 聚类
        s1 = time.time()
        cluster_pair = read_cluster_pair(cluster_file)
        cluster_pair = set(cluster_pair)
        cluster_dict = center_cluster(cluster_pair)
        clusters = [[c] + list(set(cluster_dict[c])) for c in cluster_dict if len(cluster_dict[c]) >= 0]
        s2 = time.time()
        print("Sequences clustering successfully!")
        # print('Clustering runtime: ', round(s2 - s1, 2), 's')

        # # 多序列比对
        # s3 = time.time()
        # msa_results = clusters_msa(clusters, reads, 15)
        # s4 = time.time()
        # candidates = generate_candidates(msa_results)
        # print("MSA and candidates generate successfully!")

        print('**********************************Calculate performance and show results!******************************')
        print('Loading......\n')

        print('==========Number of sequences and clusters==========')
        print("The number of original sequences is: ", len(original_sequences))
        print("The number of sequencing sequences is: ", len(reads))
        print('The number of Final Clusters is: ', len(clusters))

        print('==========Calculate runtime for each module==========')
        print('Clustering runtime: ', round(s2 - s1, 2), 's')
        # print("Multiple Sequence Alignment runtime: ", round(s4 - s3, 2), 's')

        print('==========Calculate Clustering external assessment indicators==========')
        true_labels = getting_true_labels(sequencing_file, correct_length)  # true labels
        pred_labels = getting_cluster_labels(clusters, len(reads))  # cluster labels
        # print(len(true_labels), true_labels)
        # print(len(pred_labels), pred_labels)
        nmi = compute_NMI(true_labels, pred_labels)
        print('NMI: ', round(nmi, 3))
        ami = compute_AMI(true_labels, pred_labels)
        print('AMI: ', round(ami, 3))
        ari = compute_ARI(true_labels, pred_labels)
        print('ARI: ', round(ari, 3))
        fmi = compute_FMI(true_labels, pred_labels)
        print('FMI: ', round(fmi, 3))
        homo = compute_HOMO(true_labels, pred_labels)
        print('HOMO: ', round(homo, 3))
        comp = compute_COMP(true_labels, pred_labels)
        print('COMP: ', round(comp, 3))
        v_measure = compute_V_measure(true_labels, pred_labels)
        print('V_measure: ', round(v_measure, 3))
        purity = compute_purity(true_labels, pred_labels)
        print('Purity: ', round(purity, 3))
        accuracy = compute_accuracy(true_labels, pred_labels)
        print('Accuracy: ', round(accuracy, 3))

        # print('==========Calculate fraction of recovered sequences==========')
        # fraction, redundancy = fraction_recovered(candidates, original_sequences)
        # print("The fraction of recovered sequences: ", round(fraction, 4))
        # print("The redundancy of recovered sequences is: ", round(redundancy, 3))
        #
        # print('==========Calculate sequences reconstruction rate==========')
        # re_rate = reconstruction_rate(candidates, original_sequences)
        # print('The average sequences reconstruction rate is: ', round(re_rate, 4))

        # print('==========Calculate Clustering internal assessment indicators==========')
        # intra_distances = compute_intra_cluster_levenshtein(clusters, reads)
        # intra_distances = sum(intra_distances) / len(intra_distances)
        # print('The Intra-cluster average Levenshtein Distance is: ', intra_distances)
        # inter_distances = compute_inter_cluster_levenshtein(clusters, reads)
        # inter_distances = sum(inter_distances) / len(inter_distances)
        # print('The Inter-cluster average Levenshtein Distance is: ', inter_distances)
    #     total_fraction.append(fraction)
    #     total_re_rate.append(re_rate)
    # print('Average fraction of recovered sequences: ', round(sum(total_fraction) / len(total_fraction), 4))
    # print('Average sequences reconstruction rate: ', round(sum(total_re_rate) / len(total_re_rate), 4))
