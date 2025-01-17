from sequencing_preprocessing import *
from lsh_sketch import *
from greedy_clustering import *
from cluster_merging_refinement import *
from clusters_msa import *
from compute_performances import *
from tqdm import tqdm
import time

"""
1、序列排序与预处理（根据序列计数和与原长度的偏差进行排序，对序列前后引物或高错误率区域进行修剪）
2、哈希映射（将序列的局部区域划分为若干子片段，并通过LSH将其映射为哈希值）
3、贪婪聚类（由于序列上存在indels错误，对序列进行贪婪聚类时，通过哈希差异字典和哈希草图偏移进行模糊匹配）
4、更新簇中心与簇的合并（通过计算簇中所有序列对之间的距离来更新簇中心，通过计算簇中代表序列之间的相似性将单例簇与小尺寸簇归并到好的簇中）
5、多序列比对生成候选序列（用Muscle5对每个归并后的簇进行比对，并通过投票得分来确定候选序列）
"""


def main():
    encoding_file = '../datasets/test_data/test_origin_50.fasta'  # 编码后的原始数据
    sequencing_file = '../datasets/test_data/test_cluster_50_with_labels.fasta'  # 测序后的下机数据
    # encoding_file = '../datasets/encoding/CNR_encoding.fasta'  # 编码后的原始数据
    # sequencing_file = '../datasets/sequencing/CNR_sequencing_with_labels.fasta'  # 测序后的下机数据
    k = 14  # k-mer的大小 14
    gap = 2  # 每个k-mer之间的间隔 2
    drift = 2  # 横向漂移范围 2
    sketch_size = 5  # 每条序列映射成哈希草图的尺寸 5
    star_position = 0  # 序列切片的起始位置 0
    good_bad = 5  # 区分簇的好坏的阈值,该值指代簇的尺寸 5，3
    correct_length = 110  # 序列的原始长度 100-120
    update_match_threshold = 3 / sketch_size  # 更新簇时序列间的Jaccard匹配阈值
    span = (2 * drift + k) * sketch_size + gap * (sketch_size - 1)
    end_position = star_position + span

    # encoding_file = '../datasets/encoding/I16_S2_R1_001_encoding.fasta'  # 编码后的原始数据
    # sequencing_file = '../datasets/sequencing/I16_S2_R1_001_sequencing_with_labels.fasta'  # 测序后的下机数据
    # k = 14  # k-mer的大小 14
    # gap = 2  # 每个k-mer之间的间隔 2
    # drift = 2  # 横向漂移范围 2
    # sketch_size = 2  # 每条序列映射成哈希草图的尺寸 5
    # star_position = 0  # 序列切片的起始位置 0
    # good_bad = 3  # 区分簇的好坏的阈值,该值指代簇的尺寸 3
    # correct_length = 60  # 序列的原始长度 100-120
    # update_match_threshold = 2 / sketch_size  # 更新簇时序列间的Jaccard匹配阈值
    # span = (2 * drift + k) * sketch_size + gap * (sketch_size - 1)
    # end_position = star_position + span

    # encoding_file = '../datasets/simulate_data/origin.fasta'  # 编码后的原始数据
    # sequencing_file = '../datasets/simulate_data/error_0.05/simulate_error0.05_iteration1.fasta'  # 测序后的下机数据
    # encoding_file = '../datasets/encoding/P10_5_BDDP210000009_1A_encoding.fasta'  # 编码后的原始数据
    # sequencing_file = '../datasets/sequencing/P10_5_BDDP210000009_1A_sequencing_with_labels.fasta'  # 测序后的下机数据
    # k = 15    # k=15  nmi_1=0.979  num=238305+121965=360270
    # gap = 7
    # drift = 2
    # sketch_size = 5
    # star_position = 40  # [34:-18]

    # k = 14  # k=14  nmi_1=0.953  num=174585+35012=209597
    # gap = 5
    # drift = 2
    # sketch_size = 6
    # star_position = 39  # [34:-18]

    # k = 16  # k=16  nmi_1=
    # gap = 4
    # drift = 2
    # sketch_size = 6
    # star_position = 34  # [34:-18]

    # good_bad = 3  # 区分簇的好坏的阈值,该值指代簇的尺寸 3
    # correct_length = 200  # 序列的原始长度 34-182
    # update_match_threshold = 3 / sketch_size  # 更新簇时序列间的Jaccard匹配阈值
    # span = (2 * drift + k) * sketch_size + gap * (sketch_size - 1)
    # end_position = star_position + span

    print('****************************************Starting the main programme!***************************************')
    iteration = 1
    total_fraction = []
    total_re_rate = []
    result_file = 'result.txt'
    output_file = '../results/clustering_result.txt'
    output_file_new = '../results/clustering_result_new.txt'
    for it in range(iteration):
        print('1.Sequencing preprocessing AND Read sequencing files===================================================')
        s1 = time.time()
        sequencing_sorted_file = sequencing_preprocessing(sequencing_file, correct_length)
        # sequencing_sorted_file = '../datasets/sequencing/P10_5_BDDP210000009_1A_sequencing_with_labels_sorted.fasta'
        s2 = time.time()  # 序列预处理（排序）时间
        original_sequences, reads = read_file(encoding_file, sequencing_sorted_file, correct_length)
        print("Data read successfully!")
        true_labels = getting_true_labels(sequencing_sorted_file, correct_length)  # true labels

        print('2.Generate hash signatures for each sequence===========================================================')
        s3 = time.time()
        # reads_lsh_index = [segment_lsh_index(seq[star_position:end_position], k, gap, drift) for seq in reads]
        reads_lsh_sketches = []
        for seq in tqdm(reads, desc='Generate hash sketches', ncols=100):
            read_lsh_index = segment_lsh_index(seq[star_position:end_position], k, gap, drift)
            reads_lsh_sketches.append(read_lsh_index)
        s4 = time.time()  # 序列草图构建时间
        diff_list = generate_difference_list(k)
        # reads_lsh_sketches = list_to_tuple(reads_lsh_index)
        print("Sketches build successfully!")

        print('3.Greedy clustering====================================================================================')
        s5 = time.time()
        clusters_dict = clustering_variable(reads_lsh_sketches, diff_list, sketch_size, drift)
        s6 = time.time()  # 贪婪聚类时间
        print("Sequences clustering successfully!")

        print('4.Update clusters representative and Cluster merging and refinement====================================')
        good_clusters = [[c] + list(set(clusters_dict[c])) for c in clusters_dict if len(clusters_dict[c]) >= good_bad]
        bad_clusters = [[c] + list(set(clusters_dict[c])) for c in clusters_dict if len(clusters_dict[c]) < good_bad]

        ##########
        total_clusters = good_clusters + bad_clusters
        # print(f"总簇的数量：{len(total_clusters)}\n好簇的数量：{len(good_clusters)}\n坏簇的数量：{len(bad_clusters)}")
        pred_labels_1 = getting_cluster_labels(total_clusters, len(reads))  # cluster labels
        nmi_1 = compute_NMI(true_labels, pred_labels_1)
        # print('NMI: ', round(nmi_1, 3))
        with open(result_file, 'w') as rfile:
            rfile.write(f'Total clusters: {len(total_clusters)}\n')
            rfile.write(f'Good clusters: {len(good_clusters)}\n')
            rfile.write(f'Bad clusters: {len(bad_clusters)}\n')
            rfile.write(f'Init clustering NMI: {nmi_1}\n')
        ##########

        s7 = time.time()
        for cluster in tqdm(good_clusters, desc='Update representation', ncols=100):
            cluster = update_representation(cluster, reads_lsh_sketches, diff_list, update_match_threshold)
        s8 = time.time()  # 更新簇中心时间
        merge_clusters = cluster_merging(good_clusters, reads_lsh_sketches, diff_list, sketch_size, drift)
        s9 = time.time()  # 簇的合并时间
        # final_clusters = good_clusters + bad_clusters  # 合并和细化都没有
        # final_clusters = merge_clusters + bad_clusters  # 只有合并
        # final_clusters = cluster_refinement(good_clusters, bad_clusters, reads_lsh_sketches, diff_list, sketch_size, drift)  # 只有细化
        final_clusters = cluster_refinement(merge_clusters, bad_clusters, reads_lsh_sketches, diff_list, sketch_size,
                                            drift)
        s10 = time.time()  # 簇的细化时间
        print("Cluster update representative and merging successfully!")

        print(final_clusters)
        # write_fasta(final_clusters, reads, output_file)
        write_fasta_new(final_clusters, reads, sequencing_sorted_file, output_file_new)

        print('==========Number of sequences and clusters==========')
        print("The number of original sequences is: ", len(original_sequences))
        print("The number of sequencing sequences is: ", len(reads))
        print('The number of Good Clusters is: ', len(good_clusters))
        print('The number of Bad Clusters is: ', len(bad_clusters))
        print('The number of merged Clusters is: ', len(merge_clusters))
        print('The number of Final Clusters is: ', len(final_clusters))

        print('==========Calculate runtime for each module==========')
        print('Sequences preprocessing runtime: ', round(s2 - s1, 2), 's')
        print('Constructing LSH Sketches runtime: ', round(s4 - s3, 2), 's')
        print('Initial Clustering runtime: ', round(s6 - s5, 2), 's')
        print("Update clusters representative runtime: ", round(s8 - s7, 2), 's')
        print("Cluster merging runtime: ", round(s9 - s8, 2), 's')
        print("Cluster refinement runtime: ", round(s10 - s9, 2), 's')
        print('Clustering total runtime: ', round((s6 - s5) + (s8 - s7) + (s9 - s8) + (s10 - s9), 2), 's')

        print('==========Calculate Clustering external assessment indicators==========')
        # true_labels = getting_true_labels(sequencing_sorted_file, correct_length)  # true labels
        pred_labels = getting_cluster_labels(final_clusters, len(reads))  # cluster labels
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
        # v_measure = compute_V_measure(true_labels, pred_labels)
        # print('V_measure: ', round(v_measure, 3))
        # purity = compute_purity(true_labels, pred_labels)
        # print('Purity: ', round(purity, 3))
        # accuracy = compute_accuracy(true_labels, pred_labels)
        # print('Accuracy: ', round(accuracy, 3))

        # print('==========Calculate Clustering internal assessment indicators==========')
        # intra_distances = compute_intra_cluster_levenshtein(final_clusters, reads)
        # intra_distances = sum(intra_distances) / len(intra_distances)
        # print('The Intra-cluster average Levenshtein Distance is: ', intra_distances)
        # inter_distances = compute_inter_cluster_levenshtein(final_clusters, reads)
        # inter_distances = sum(inter_distances) / len(inter_distances)
        # print('The Inter-cluster average Levenshtein Distance is: ', inter_distances)

        with open(result_file, 'a') as rfile:
            rfile.write(f'The number of original sequences is: {len(original_sequences)}\n')
            rfile.write(f'The number of sequencing sequences is: {len(reads)}\n')
            rfile.write(f'The number of Good Clusters is: {len(good_clusters)}\n')
            rfile.write(f'The number of Bad Clusters is: {len(bad_clusters)}\n')
            rfile.write(f'The number of merged Clusters is: {len(merge_clusters)}\n')
            rfile.write(f'The number of Final Clusters is: {len(final_clusters)}\n')
            rfile.write(f'Clustering time: {round((s6 - s5) + (s8 - s7) + (s9 - s8) + (s10 - s9), 2)} s \n')
            rfile.write(f'NMI: {round(nmi, 3)}\n')
            rfile.write(f'AMI: {round(ami, 3)}\n')
            rfile.write(f'ARI: {round(ari, 3)}\n')
            rfile.write(f'FMI: {round(fmi, 3)}\n')
            rfile.write(f'HOMO: {round(homo, 3)}\n')
            rfile.write(f'COMP: {round(comp, 3)}\n')

        # print('5.Multiple Sequence Alignment AND Majority voting for Candidates=======================================')
        # s11 = time.time()
        # msa_results = clusters_msa(final_clusters, reads, 15)
        # s12 = time.time()  # 多序列比对时间
        # candidates = generate_candidates(msa_results)
        # print("MSA and candidates generate successfully!")
        # print("Multiple Sequence Alignment runtime: ", round(s12 - s11, 2), 's')
        # with open(result_file, 'a') as rfile:
        #     rfile.write(f'MSA time: {round(s12 - s11, 2)} s \n')
        #
        # print('==========Calculate fraction of recovered sequences==========')
        # fraction, redundancy = fraction_recovered(candidates, original_sequences)
        # print("The fraction of recovered sequences: ", round(fraction, 4))
        # print("The redundancy of recovered sequences is: ", round(redundancy, 3))
        # with open(result_file, 'a') as rfile:
        #     rfile.write(f'Fraction recovered: {round(fraction, 4)}\n')
        #
        # print('==========Calculate sequences reconstruction rate==========')
        # re_rate = reconstruction_rate(candidates, original_sequences)
        # print('The sequences reconstruction rate is: ', round(re_rate, 4))
        # with open(result_file, 'a') as rfile:
        #     rfile.write(f'Reconstruction rate: {round(re_rate, 4)}\n')
        #
        # ami = compute_AMI(true_labels, pred_labels)
        # print('AMI: ', round(ami, 3))

        #     total_fraction.append(fraction)
        #     total_re_rate.append(re_rate)
        # print('Average fraction of recovered sequences: ', round(sum(total_fraction) / len(total_fraction), 4))
        # print('Average sequences reconstruction rate: ', round(sum(total_re_rate) / len(total_re_rate), 4))
        # print(final_clusters)

        true_labels_new = getting_true_labels_new(sequencing_file)
        pred_labels_new = getting_cluster_labels_new(output_file_new)
        # print(len(true_labels_new), true_labels_new)
        # print(len(pred_labels_new), pred_labels_new)
        nmi_new = compute_NMI(true_labels_new, pred_labels_new)
        print('NMI_new: ', round(nmi_new, 3))
        ami_new = compute_AMI(true_labels_new, pred_labels_new)
        print('AMI_new: ', round(ami_new, 3))
        ari_new = compute_ARI(true_labels_new, pred_labels_new)
        print('ARI_new: ', round(ari_new, 3))
        fmi_new = compute_FMI(true_labels_new, pred_labels_new)
        print('FMI_new: ', round(fmi_new, 3))
        homo_new = compute_HOMO(true_labels_new, pred_labels_new)
        print('HOMO_new: ', round(homo_new, 3))
        comp_new = compute_COMP(true_labels_new, pred_labels_new)
        print('COMP_new: ', round(comp_new, 3))


if __name__ == '__main__':
    main()
