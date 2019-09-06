from bio_calculator import Calculator
from pattern_dna import PatternDNA

read_file_path = "./data_set/X174.0.00.protein.fa"


def get_source_pool(path):
    source_pool = []
    file = open(path, "r")
    lines = file.read().splitlines()
    for index in range(int(len(lines) / 2)):
        source_pool.append(PatternDNA(indices=[str(lines[index * 2][1:])], protein=lines[index * 2 + 1]))

    return source_pool


def get_overlap_matrix(source_pool):
    total_calculate_count = int((len(source_pool) * (len(source_pool) - 1)) / 2)
    overlap_matrix = [[0 for col in range(len(source_pool))] for row in range(len(source_pool) - 1)]
    calculator_map = {}
    greedy = []
    current_count = 0
    for row in range(len(overlap_matrix)):
        for col in range(row + 1, len(overlap_matrix[row])):
            calculator = Calculator()
            calculator.calculate_overlap(source_pool[row], source_pool[col])
            overlap_length = calculator.get_overlap_length()

            if overlap_length > 0:
                greedy.append([overlap_length, row, col])
                overlap_matrix[row][col] = overlap_length
                calculator_map[str([row, col])] = calculator

            current_count += 1
            print("\rDetect = " + str(current_count) + " | Total = " + str(total_calculate_count), end=" ")

    print("\nresult = ")
    for row in range(len(overlap_matrix)):
        a = str(overlap_matrix[row])[1:-1]
        print(a)


if __name__ == "__main__":
    pool = get_source_pool("./data_set/X174.0.00.protein.fa")
    get_overlap_matrix(pool)