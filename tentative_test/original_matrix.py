import numpy

from utils import accessor

if __name__ == "__main__":
    indices, data_domains, whole_dna_seq = accessor.read_information_by_gb("../data_set/Ecoli_K-1_MG1655.gb")
    overlap_matrix = [[0 for col in range(len(indices))] for row in range(len(indices) - 1)]

    current = 0
    total = int((len(indices) * (len(indices) - 1)) / 2)
    for row in range(len(indices)):
        for col in range(row + 1, len(indices)):
            dna_1 = numpy.array([index for index in range(data_domains[row][1], data_domains[row][2])])
            dna_2 = numpy.array([index for index in range(data_domains[col][1], data_domains[col][2])])
            overlap = len(numpy.intersect1d(dna_1, dna_2))
            overlap_matrix[row][col] = overlap
            current += 1
            print("\r" + "Detect: " + str(current) + " (" + str(total) + ")" +
                  " || (" + str(row) + ", " + str(col) + ") = " + str(overlap), end=" ")

    # accessor.format_output("original matrix", overlap_matrix)
    with open("../output/Ecoli_K-1_MG1655.csv", "w", encoding="utf-8") as save_file:
        for row in range(len(overlap_matrix)):
            row_data = str(overlap_matrix[row])[1: -1]
            save_file.write(row_data + "\n")

