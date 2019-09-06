"""
Name: PatternDNA

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1
"""

from bio_calculator import Calculator

from pattern_dna import PatternDNA

read_file_path = "./data_set/X174.0.00.protein.fa"

if __name__ == "__main__":
    # init process
    source_pool = []
    target_pool = []

    file = open(read_file_path, "r")
    lines = file.read().splitlines()
    for index in range(int(len(lines) / 2)):
        source_pool.append(PatternDNA(indices=[str(lines[index * 2][1:])], protein=lines[index * 2 + 1]))

    print("Source pool has " + str(len(source_pool)) + " Pattern DNA.")

    # calculating process
    calculating_round = 1
    while len(source_pool) > 1:
        print("Round " + str(calculating_round))

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
                    calculator_map[[row, col]] = calculator

                current_count += 1
                print("\r" + "Detect = " + str(current_count) + " | Total = " + str(total_calculate_count), end=" ")

            # no overlap
            if sum(overlap_matrix[row]) == 0:
                target_pool.append(source_pool[row])

        print("result = ")
        print(overlap_matrix)

        greedy.sort(reverse=True)

        available = [True for index in range(len(source_pool))]
        source_pool = []

        for fuse_information in greedy:
            index_1 = fuse_information[1]
            index_2 = fuse_information[2]
            if available[index_1] and available[index_2]:
                source_pool.append(calculator_map.get([index_1, index_2]).get_fuse_dna())
                available[index_1] = False
                available[index_2] = False

        print("Current source pool has " + str(len(source_pool)) + " Pattern DNA.")
        calculating_round += 1

    for dna in source_pool:
        target_pool.append(dna)

    # outputting process
    for dna in target_pool:
        print(dna)
