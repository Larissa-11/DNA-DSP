import json

# def read_dna_file(path):
#     """
#     introduction: Reading DNA sequence set from documents.
#
#     :param path: File path.
#                   Type: string
#
#     :return dna_sequences: A corresponding DNA sequence string in which each row acts as a sequence.
#                            Type: one-dimensional list(string)
#
#     """
#     dna_sequences = []
#     with open(path, "r") as file:
#        # Read current file by line
#         lines = file.readlines()
#         for index in range(len(lines)):
#             line = lines[index]
#             dna_sequences.append([line[col] for col in range(len(line) - 1)])
#
#     return dna_sequences


def read_dna_file(path):
    """
    introduction: Reading DNA sequence set from documents.

    :param path: File path.
                  Type: string

    :return dna_sequences: A corresponding DNA sequence string in which each row acts as a sequence.
                           Type: one-dimensional list(string)
    """
    dna_sequences = []
    with open(path, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()  # 去除行首尾的空白符
            if line:  # 检查是否为空白行
                dna_sequences.append(list(line))
    return dna_sequences


def write_dna_file(path, dna_sequences):
    """
    introduction: Writing DNA sequence set to documents.

    :param path: File path.
                  Type: string

    :param dna_sequences: Generated DNA sequences.
                          Type: one-dimensional list(string)
    """

    with open(path, "w") as file:
        for row in range(len(dna_sequences)):
            file.write("".join(dna_sequences[row]) + "\n")
    return dna_sequences


def load_dict_from_file(filepath):
    with open(filepath, 'r') as file:
        dictionary = json.load(file)
    return dictionary
