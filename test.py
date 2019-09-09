from common.fuser import Fuser
from common.dna import PatternDNA


if __name__ == "__main__":
    source_pool = []
    file = open("./data_set/X174.0.00.protein.fa", "r")
    lines = file.read().splitlines()
    for index in range(int(len(lines) / 2)):
        source_pool.append(PatternDNA(indices=[str(lines[index * 2][1:])], protein=lines[index * 2 + 1]))

    total_calculate_count = int((len(source_pool) * (len(source_pool) - 1)) / 2)
    overlap_matrix = [[0 for col in range(len(source_pool))] for row in range(len(source_pool) - 1)]
    current_count = 0
    for row in range(len(source_pool)):
        for col in range(row + 1, len(source_pool)):
            fuser = Fuser()
            fuser.calculate_overlap(source_pool[row], source_pool[col])
            overlap_length = fuser.get_overlap_length()

            current_count += 1

            if overlap_length > 0:
                overlap_matrix[row][col] = overlap_length
                print(str(row) + ", " + str(col) + ", " + str(overlap_length) +
                      ", " + str(current_count) + ", " + str(total_calculate_count))


    print(overlap_matrix)
