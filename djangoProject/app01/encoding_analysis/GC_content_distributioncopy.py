import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np

from app01.encoding_analysis import data_handle


def compute_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = gc_count / len(sequence) * 100
    return gc_content


# 生成 DNA 序列并存储到不同文件中
encoding_methods = ['./DNAFountain.dna','./DNA-QLC.dna', './0801.dna']
gc_contents_all = []
for encoding_method in encoding_methods:
    dna_sequences = data_handle.read_dna_file(encoding_method)
    gc_contents = [compute_gc_content(sequence) for sequence in dna_sequences]
    gc_contents_all.append(gc_contents)

# 设置样式
sns.set_style("white")

# 绘制直方图
fig, axs = plt.subplots(nrows=3, figsize=(8, 8), sharex=True, sharey=True)
plt.subplots_adjust(top=0.95, bottom=0.1, left=0.12, right=0.95, hspace=0.1)

for idx, ax in enumerate(axs):
    sns.histplot(gc_contents_all[idx], kde=False, ax=ax,stat='percent')
    ax.set_xlim(40, 60)
    ax.set_xticks(np.arange(40, 61, 5))
    ax.text(0.05, 0.9, encoding_methods[idx].split('/')[-1].split('.')[0], transform=ax.transAxes, ha='left', va='top',
            fontsize='12')

plt.xlabel('GC Content (%)')
# 保存图形
plt.savefig('GC_content.png')

# 选择随机 DNA 片段并计算不同长度的 GC 含量
sequence_lengths = list(range(6, 151, 6))

gc_data = {'Sequence Length': [], 'GC Content (%)': [], 'Encoding Method': []}

for encoding_method in encoding_methods:
    dna_sequence = data_handle.read_dna_file(encoding_method)
    random_sequence = random.choice(dna_sequence)

    for length in sequence_lengths:
        sequence_subset = random_sequence[:length]
        gc_content = compute_gc_content(sequence_subset)
        gc_data['Sequence Length'].append(length)
        gc_data['GC Content (%)'].append(gc_content)
        gc_data['Encoding Method'].append(encoding_method.split('/')[-1].split('.')[0])

gc_df = pd.DataFrame(gc_data)

plt.figure(figsize=(6, 6))
encoding_colors = {
    # 将0801编码方式对应的线条颜色设置为红色
    'DNAFountain': 'orange',
    'DNA-QLC': 'green',
    '0801': '#1f77b4'
}
# 绘制线图，并根据编码方式指定线条颜色
sns.lineplot(data=gc_df, x='Sequence Length', y='GC Content (%)', hue='Encoding Method', palette=encoding_colors,
             errorbar=None)
plt.xlabel("Sequence Length (nt)")
plt.ylabel("GC Content (%)")
plt.legend(loc='upper right')
plt.grid(True)
plt.subplots_adjust(top=0.95)
plt.show()
plt.savefig('Local_GC_content.png')
