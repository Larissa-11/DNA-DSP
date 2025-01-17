from nupack import *
import os

DFT_CONCENTRATION = 1e-1
MODEL = Model(material='dna', celsius=37, sodium=1.0, magnesium=0.0, ensemble='stacking')
input_files = ['./1.dna', './2.dna', './3.dna']
mfe=[]
for file in input_files:
    score_all = []
    dna_sequences = []
    with open(file, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()  # 去除行首尾的空白符
            if line:  # 检查是否为空白行
                dna_sequences.append(list(line))
    for i in dna_sequences:
        seq = "".join(i)
        seq1 = Strand(seq, name="a")
        my_complex = Complex([seq1], name="b")
        tube1 = Tube({seq1: DFT_CONCENTRATION}, complexes=SetSpec(max_size=1, include=[my_complex]), name="tube1")
        single_results = tube_analysis([tube1], model=MODEL, compute=['pfunc', 'pairs', 'mfe', 'sample', 'subopt'],
                                       options={'num_sample': 2, 'energy_gap': 0.5})
        score = single_results[my_complex].mfe[0].energy
        score_all.append(score)
    mfe.append(score_all)
