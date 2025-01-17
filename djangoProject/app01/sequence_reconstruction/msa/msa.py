import operator
import subprocess
from Bio import AlignIO
from tqdm import tqdm
import os


def process_byte_stream(byte_stream):
    sequences = []
    clusters = []
    current_cluster = []

    # 将字节流按行进行处理
    lines = byte_stream.split(b'\n')
    for line in lines:
        line = line.strip()
        if line.startswith(b'>'):
            # 遇到新的 Cluster 标识，将当前的 cluster 加入 clusters
            if current_cluster:
                clusters.append(current_cluster)
                current_cluster = []
        elif line and not line.startswith(b'='):
            # 处理 ATCG 开头的行，构建 sequences
            if line.startswith(b'A') or line.startswith(b'T') or line.startswith(b'C') or line.startswith(b'G'):
                sequences.append(line.decode())
                current_cluster.append(len(sequences) - 1)

    # 添加最后一个 cluster
    if current_cluster:
        clusters.append(current_cluster)

    return sequences, clusters


def format_sequences_as_fasta(sequences):
    formatted_fasta = []
    for i, seq in enumerate(sequences):
        # seq_id = f">Seq_{i + 1}"
        formatted_fasta.append(f"{seq}")
    return '\n'.join(formatted_fasta)


def msa(cluster):
    # script_dir = os.path.dirname(os.path.abspath(__file__))  # 获取当前脚本所在目录的绝对路径
    # input_alignment = os.path.join(script_dir, "clm.fasta")  # 创建输入文件的完整路径
    # output_alignment = os.path.join(script_dir, "clmout.fasta")
    input_alignment = "app01/sequence_reconstruction/msa/clm.fasta"
    output_alignment = "app01/sequence_reconstruction/msa/clmout.fasta"
    with open(input_alignment, "w") as file:
        for i, c in enumerate(cluster):
            file.write(f">S{i}\n")
            file.write(c)
            file.write("\n")
    with open(output_alignment, "w") as file:
        pass
    # 构建执行多序列比对工具的命令（可修改）
    # msa_exe = os.path.join(script_dir, 'muscle5.exe')
    msa_exe = "app01/sequence_reconstruction/msa/muscle5.1.linux_intel64"
    msa_command = [msa_exe, "-align", input_alignment, "-output", output_alignment]

    try:
        subprocess.run(msa_command, check=True)
        # print("Alignment completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running: {e}")
    msa_sequences = AlignIO.read(output_alignment, "fasta")
    aligned_cluster = []
    for i in msa_sequences:
        aligned_cluster += [i.seq]
    return aligned_cluster


def clusters_msa(clusters, sequences):
    results = []
    for i, cluster_indexes in enumerate(tqdm(clusters, ncols=100)):
        cluster = [sequences[index] for index in cluster_indexes]
        aligned_cluster = msa(cluster)
        results.append(aligned_cluster)
        # if i % 1000 == 0:
        #     print("%", round(i * 100 / len(clusters), 2), "of the clusters are aligned.")
    return results


def sequence_voting(m, weight):
    res = ""
    for i in range(len(m[0])):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0, 'N': 0}
        for j in range(len(m)):
            counts[m[j][i]] += 1
        counts['-'] *= weight
        mv = max(counts.items(), key=operator.itemgetter(1))[0]
        if mv != '-':
            res += mv
    return res


def generate_candidates(results):
    candidates = []
    for m in results:
        candidates.append(sequence_voting(m, weight=0.5))
    return candidates




