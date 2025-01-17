from Chamaeleo.utils.monitor import Monitor
import numpy as np


def read_bits_from_byte(file_byte, segment_length=120, need_logs=False):
    monitor = Monitor()
    if need_logs:
        print("Read binary matrix from file: " + file_byte)

    matrix, values = [], np.frombuffer(file_byte, dtype=np.uint8)
    for current, value in enumerate(values):
        matrix += list(map(int, list(str(bin(value))[2:].zfill(8))))
        if need_logs:
            monitor.output(current + 1, len(values))
    if len(matrix) % segment_length != 0:
        matrix += [0] * (segment_length - len(matrix) % segment_length)

    matrix = np.array(matrix)
    matrix = matrix.reshape(int(len(matrix) / segment_length), segment_length)

    if need_logs:
        print("There are " + str(len(values) * 8) + " bits in the inputted file. "
              + "Please keep this information in mind if you do not consider storing the model in serialization!")

    return matrix.tolist(), len(values) * 8


def read_dna_byte(file_byte, need_logs=False):
    monitor = Monitor()
    dna_sequences = []

    content_str = file_byte.decode('utf-8')
    content_list = content_str.split("\n")
    for index, line in enumerate(content_list):
        dna_sequences.append(list(line.replace("\n", "")))

        if need_logs:
            monitor.output(index + 1, len(content_list))

    return dna_sequences
