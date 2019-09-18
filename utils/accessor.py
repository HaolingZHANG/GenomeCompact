import math
from Bio import SeqIO


def read_information_by_fa(path):
    indices = []
    data_domains = []

    with open(path, "r") as fasta_file:
        fasta_data = fasta_file.read().splitlines()

        for index in range(int(len(fasta_data) / 2)):
            indices.append(str(fasta_data[index * 2][1:]))
            data_domains.append(list(fasta_data[index * 2 + 1]))

    return indices, data_domains


def read_information_by_gb(path, data_type="CDS"):
    indices = []
    data_domains = []
    with open(path, "r") as gb_file:
        genebank_data = SeqIO.read(gb_file, "genbank")
        whole_dna_seq = list(genebank_data.seq)

        for feature in genebank_data.features:
            if feature.type == data_type:
                indices.append(feature.qualifiers.get('label'))
                data_domains.append([
                    feature.location.strand,
                    feature.location.start,
                    feature.location.end
                ])

    return indices, data_domains, whole_dna_seq


def write_information_to_fa(path, indices, data_domains):
    pass


def write_information_to_gb(path):
    pass


def format_output(title, matrix):
    max_value = -1
    for row_data in matrix:
        if max_value < max(row_data):
            max_value = max(row_data)
    message = ""
    if max_value < 10:
        for row in range(len(matrix)):
            for col in range(len(matrix[row])):
                data = str(matrix[row][col])
                message += " " + data
            message += "\n"
    else:
        format_length = math.ceil(math.log(max_value, 10)) + 2

        for row in range(len(matrix)):
            for col in range(len(matrix[row])):
                data = str(matrix[row][col])
                message += (" " * (format_length - len(data))) + data
            message += "\n"

    print(title + ": ")
    print(message)
