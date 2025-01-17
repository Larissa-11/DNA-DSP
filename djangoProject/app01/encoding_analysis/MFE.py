from nupack import *

class Mfe:
    def __init__(self):
        self.DFT_CONCENTRATION = 1e-1
        self.MODEL = Model(material='dna', celsius=37, sodium=1.0, magnesium=0.0, ensemble='stacking')
    def mfe_calculation(self, dna_sequences_list):
        mfe = []
        for dna_sequences in dna_sequences_list:
            print("9++执行***")
            score_all = []
            for i in dna_sequences:
                # print("9++++++++++++++++++++++++++执行***")
                seq = "".join(i)
                # print(seq)
                seq1 = Strand(seq, name="a")
                # print(seq1)
                my_complex = Complex([seq1], name="b")
                # print(my_complex)
                tube1 = Tube({seq1: self.DFT_CONCENTRATION}, complexes=SetSpec(max_size=1, include=[my_complex]), name="tube1")
                # print(tube1)
                single_results = tube_analysis([tube1], model=self.MODEL, compute=['pfunc', 'pairs', 'mfe', 'sample', 'subopt'],
                                               options={'num_sample': 2, 'energy_gap': 0.5})
                # print(single_results)
                score = single_results[my_complex].mfe[0].energy
                # print(score)
                score_all.append(score)
            mfe.append(score_all)
        return mfe
