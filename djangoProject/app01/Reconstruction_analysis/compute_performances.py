import re
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn import metrics
import Levenshtein
from tqdm import tqdm


# 指标1：聚类内部评估指标
def levenshtein_distance(seq1, seq2):
    return Levenshtein.distance(seq1, seq2)


def compute_intra_cluster_levenshtein(clusters, reads):  # 簇内平均Levenshtein距离
    average_distances = []
    for cluster in clusters:
        distance = []
        for i in range(len(cluster)):
            for j in range(i + 1, len(cluster)):
                distance.append(levenshtein_distance(reads[cluster[i]], reads[cluster[j]]))
        if len(distance) != 0:
            average_distance = sum(distance) / len(distance)
        else:
            average_distance = 0
        average_distances.append(average_distance)
    return average_distances


def compute_inter_cluster_levenshtein(clusters, reads):  # 簇间平均Levenshtein距离
    average_distances = []
    for i in range(len(clusters)):
        distance = 0
        for j in range(len(clusters)):
            if i != j:
                distance += levenshtein_distance(reads[clusters[i][0]], reads[clusters[j][0]])
        average_distance = distance / len(clusters) - 1
        average_distances.append(average_distance)
    return average_distances


def compute_SI(x, labels):  # 轮廓系数（Silhouette Index）
    si = metrics.silhouette_score(x, labels)
    return si


def compute_CHI(x, labels):  # Calinski-Harabaz指数（Calinski-Harabaz Index）
    chi = metrics.calinski_harabasz_score(x, labels)
    return chi


def compute_NCC(x, labels):  # Internal and external validation measures (NCC)
    # NCC 测量的是簇内距离和簇间距离的组合。簇内距离是指一个簇内物体之间的距离。作者将两个对象的 “群内协议”（Sintra）定义为两个对象之间的属性数量和群内距离之差。
    m = x.shape[0]
    n = x.shape[1]
    Y = np.zeros((m, m))
    for r in range(m):
        for s in range(m):
            if labels[r] == labels[s]:
                Y[r, s] = 1
    drs = np.zeros((m, m))
    for r in range(m):
        for s in range(m):
            for att in range(n):
                if x[r, att] != x[s, att]:
                    drs[r, s] += 1
    ncc = 0.0
    for r in range(m):
        for s in range(m):
            if r != s:
                ncc += (n - 2 * drs[r, s]) * Y[r, s] + drs[r, s]
    return ncc


# 指标2：聚类外部评估指标
def getting_cluster_labels(clusters, length):  # 获取聚类标签
    pred_dict = {}
    for i, clu in enumerate(clusters):
        for cl in clu:
            pred_dict[cl] = i
    # print(pred_dict)
    # print(len(pred_dict))
    pred_labels = []
    i = 0
    j = -1
    while i <= length - 1:
        # print(i)
        if i in pred_dict:
            pred_labels.append(pred_dict[i])
            i += 1
        else:
            pred_labels.append(j)
            j -= 1
            i += 1
    # print(pred)
    # print(len(pred))
    return pred_labels


def getting_true_labels(input_file, correct_length):  # 获取真实标签
    true_labels = []
    min_length = correct_length - 10
    max_length = correct_length + 10
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                cluster_label = int(re.search(r'Cluster(\d+)', line).group(1))
                true_labels.append(cluster_label)
            else:
                if min_length <= len(line) <= max_length:
                    continue
                else:
                    true_labels.pop()
                    continue
    return true_labels


# def getting_true_labels_new(input_file):
#     true_dict = {}
#     true_labels = []
#     i = 0
#     j = -1
#     with open(input_file, 'r') as file:
#         for line in file:
#             line = line.strip()
#             if line.startswith('>'):
#                 seq_label = int(re.search(r'Seq(\d+)', line).group(1))
#                 cluster_label = int(re.search(r'Cluster(\d+)', line).group(1))
#                 true_dict[seq_label] = cluster_label
#     while i <= len(true_dict) - 1:
#         # print(i)
#         if i in true_dict:
#             true_labels.append(true_dict[i])
#             i += 1
#         else:
#             true_labels.append(j)
#             j -= 1
#             i += 1
#     return true_labels
import re


def getting_true_labels_new(input_file):
    true_dict = {}
    true_labels = []
    i = 0
    j = -1

    # 读取 InMemoryUploadedFile 的内容
    file_content = input_file.read().decode('utf-8')
    for line in file_content.splitlines():
        line = line.strip()
        if line.startswith('>'):
            seq_label = int(re.search(r'Seq(\d+)', line).group(1))
            cluster_label = int(re.search(r'Cluster(\d+)', line).group(1))
            true_dict[seq_label] = cluster_label

    while i <= len(true_dict) - 1:
        if i in true_dict:
            true_labels.append(true_dict[i])
            i += 1
        else:
            true_labels.append(j)
            j -= 1
            i += 1

    return true_labels


# def getting_cluster_labels_new(input_file):
#     pred_dict = {}
#     pred_labels = []
#     i = 0
#     j = -1
#     with open(input_file, 'r') as file:
#         for line in file:
#             line = line.strip()
#             if line.startswith('>'):
#                 seq_label = int(re.search(r'Seq(\d+)', line).group(1))
#                 cluster_label = int(re.search(r'Cluster(\d+)', line).group(1))
#                 pred_dict[seq_label] = cluster_label
#     while i <= len(pred_dict) - 1:
#         # print(i)
#         if i in pred_dict:
#             pred_labels.append(pred_dict[i])
#             i += 1
#         else:
#             pred_labels.append(j)
#             j -= 1
#             i += 1
#     return pred_labels
def getting_cluster_labels_new(input_file):
    pred_dict = {}
    pred_labels = []
    i = 0
    j = -1

    # 读取 InMemoryUploadedFile 的内容
    file_content = input_file.read().decode('utf-8')
    for line in file_content.splitlines():
        line = line.strip()
        if line.startswith('>'):
            seq_label = int(re.search(r'Seq(\d+)', line).group(1))
            cluster_label = int(re.search(r'Cluster(\d+)', line).group(1))
            pred_dict[seq_label] = cluster_label

    while i <= len(pred_dict) - 1:
        if i in pred_dict:
            pred_labels.append(pred_dict[i])
            i += 1
        else:
            pred_labels.append(j)
            j -= 1
            i += 1

    return pred_labels

def compute_NMI(true_labels, pred_labels):  # 标准化互信息（NMI, Normalized Mutual Information）互信息是用来衡量两个数据分布的吻合程度 [0,1]
    nmi = metrics.normalized_mutual_info_score(true_labels, pred_labels)
    return nmi


def compute_AMI(true_labels, pred_labels):  # 调整互信息（AMI, Adjusted mutual information） [-1,1]
    ami = metrics.adjusted_mutual_info_score(true_labels, pred_labels)
    return ami


def compute_ARI(true_labels, pred_labels):  # 调整兰德系数（ARI, Adjusted Rand index）衡量的是两个数据分布的吻合程度 [-1,1]
    ari = metrics.adjusted_rand_score(true_labels, pred_labels)
    return ari


def compute_FMI(true_labels, pred_labels):  # FMI（Fowlkes-Mallows Index）是 Precision（精度）和 Recall（召回）的几何平均数 [0,1]
    fmi = metrics.fowlkes_mallows_score(true_labels, pred_labels)
    return fmi


def compute_HOMO(true_labels, pred_labels):  # 同质性（homogeneity）：每个群集只包含单个类的成员
    homo = metrics.homogeneity_score(true_labels, pred_labels)
    return homo


def compute_COMP(true_labels, pred_labels):  # 完整性（completeness）：给定类的所有成员都分配给同一个群集
    comp = metrics.completeness_score(true_labels, pred_labels)
    return comp


def compute_V_measure(true_labels, pred_labels):  # V-measure是同质性homogeneity和完整性completeness的调和平均数 [0,1]
    v_measure = metrics.v_measure_score(true_labels, pred_labels)
    return v_measure


def compute_purity(true_labels, pred_labels):
    clusters = np.unique(pred_labels)
    true_labels = np.reshape(true_labels, (-1, 1))
    pred_labels = np.reshape(pred_labels, (-1, 1))
    count = []

    for c in clusters:
        idx = np.where(pred_labels == c)[0]
        labels_tmp = true_labels[idx, :].reshape(-1)
        count.append(np.bincount(labels_tmp).max())

    return np.sum(count) / true_labels.shape[0]


def compute_accuracy(true_labels, pred_labels):
    unique_true_labels = np.unique(true_labels)
    unique_predicted_labels = np.unique(pred_labels)

    contingency_matrix = np.zeros((len(unique_true_labels), len(unique_predicted_labels)))

    for i, true_label in enumerate(unique_true_labels):
        true_mask = (true_labels == true_label)
        for j, predicted_label in enumerate(unique_predicted_labels):
            predicted_mask = (pred_labels == predicted_label)
            contingency_matrix[i, j] = np.sum(np.logical_and(true_mask, predicted_mask))

    row_ind, col_ind = linear_sum_assignment(-contingency_matrix)

    accuracy = contingency_matrix[row_ind, col_ind].sum() / len(true_labels)

    return accuracy


# 指标3：序列回收分数与冗余------原始序列中有多少比例成功恢复
def fraction_recovered(candidates, orig_seqs):
    count_dict = {}
    redundancy = 0

    for seq in orig_seqs:
        count_dict[seq] = 0
    for cand in candidates:
        if cand in count_dict:
            count_dict[cand] += 1

    fraction = sum([count_dict[seq] > 0 for seq in count_dict]) / len(count_dict)
    if fraction > 0:
        redundancy = sum([count_dict[seq] for seq in count_dict]) / len(count_dict) / fraction

    return fraction, redundancy


# 指标4：序列平均重建率------重建后的序列与正确的序列之间的重建相似度
def reconstruction_rate(candidates, orig_seqs):
    reconstruct_dict = {}
    for seq in orig_seqs:
        reconstruct_dict[seq] = 0
    for cand in tqdm(candidates, desc='Reconstruction rate', ncols=100):
        if cand in reconstruct_dict:
            reconstruct_dict[cand] += 1
        else:
            # 计算该序列与orig_seqs中哪条序列间的编辑距离最短，然后计算该条序列的重建率
            min_distance = float('inf')
            target_seq = None
            for seq in orig_seqs:
                distance = levenshtein_distance(cand, seq)
                if distance < min_distance:
                    min_distance = distance
                    target_seq = seq
            re_rate = 1 - min_distance / len(target_seq)
            reconstruct_dict[target_seq] += re_rate

    average_re_rate = sum([reconstruct_dict[seq] for seq in reconstruct_dict]) / len(candidates)

    return average_re_rate


if __name__ == "__main__":
    print(1)
    # input = '../datasets/test_data/test_cluster_50_with_labels_sorted.fasta'
    # true_labels = getting_true_labels(input)
    # print(true_labels)
    # print(len(true_labels))

    # true_labels = [1, 1, 1, 1, 2, 2, 2, 3, 3]
    # pred_labels = [11, 11, 4, 11, 22, 111, 22, 22, 4]
    # nmi = compute_NMI(true_labels, pred_labels)
    # ami = compute_AMI(true_labels, pred_labels)
    # ari = compute_ARI(true_labels, pred_labels)
    # fmi = compute_FMI(true_labels, pred_labels)
    # homo = compute_HOMO(true_labels, pred_labels)
    # comp = compute_COMP(true_labels, pred_labels)
    # v_measure = compute_V_measure(true_labels, pred_labels)
    # purity = compute_purity(true_labels, pred_labels)
    # accuracy = compute_accuracy(true_labels, pred_labels)
    # print('NMI:', nmi)
    # print('AMI:', ami)
    # print('ARI:', ari)
    # print('FMI:', fmi)
    # print('HOMO:', homo)
    # print('COMP:', comp)
    # print('V_measure:', v_measure)
    # print('Purity:', purity)
    # print('Accuracy:', accuracy)

    # candidates = [
    #     'GAGAAGGGGCCTCGAGTTTATCCCTAACCAAACGGTGAGCATGAGCGTTGACCCTGCATGATAGTCCGGTGGTCGAGGAGATAGACATCGGATGATTTCAACATGTAAGT',
    #     'GTGGCGACCAACTTCACTTAGATTGTCGGTCGCAGTGTTACTCGGTCGCGCGGTTGAAAGGGGTTAGCAGCACTCACGGTCCAAGCAGGCTAAGCCTTATTCCCAAACAG',
    #     'GAAGGGGCCTCGAGTTTATCCCTAACCAAACGGTGAGACATGAGCGTTGACCCTGCATGATAGTCCGGTGGTCGAGGAAATAGACATCGGAATGATTTCAACATGTAAAT',
    #     'GGCGACCAACTTCACTTAGATTGTCGATCGCAGTGTTACTCGGTCGCGCGGTTGAAAGGGGTTAGCAGCACTCACGATCCAAGCAGGCTAAGCCATATTCCCAAACAGAA']
    # orig_seqs = [
    #     'GAGAAGGGGCCTCGAGTTTATCCCTAACCAAACGGTGAGCATGAGCGTTGACCCTGCATGATAGTCCGGTGGTCGAGGAGATAGACATCGGATGATTTCAACATGTAAGT',
    #     'GTGGCGACCAACTTCACTTAGATTGTCGGTCGCAGTGTTACTCGGTCGCGCGGTTGAAAGGGGTTAGCAGCACTCACGGTCCAAGCAGGCTAAGCCTTATTCCCAAACAG']
    # average_re_rate = reconstruction_rate(candidates, orig_seqs)
    # print(average_re_rate)
