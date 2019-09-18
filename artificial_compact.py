"""
Name: PatternDNA

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1
"""

import multiprocessing

from methods.fuser import Fuser
from methods.fuser import init_indices
from methods.fuser import init_matched_group

from methods.dna import PatternDNA

read_file_path = "./data_set/X174.0.00.protein.fa"
write_file_path = "./output/pattern_dna.fa"


def get_source_pool(path):
    pool = []
    fasta_file = open(path, "r")
    proteins = fasta_file.read().splitlines()
    for protein_index in range(int(len(proteins) / 2)):
        pool.append(PatternDNA(indices=[str(proteins[protein_index * 2][1:])], protein=proteins[protein_index * 2 + 1]))

    print("Original source pool has " + str(len(pool)) + " Pattern DNAs.")

    return pool


def calculate(pool, row, col, lock, current, total):
    fuser = Fuser(init_indices(pool[row], pool[col]))
    overlap_length = fuser.calculate_overlap(init_matched_group(pool[row], pool[col]))

    lock.acquire()
    current.value += 1
    lock.release()
    print("\rdetect: " + str(current.value) + ", (" + str(total.value) +
          ") || (" + str(row) + ", " + str(col) + ") = " + str(overlap_length), end=" ")

    return [row, col, overlap_length, fuser]


# noinspection PyUnusedLocal
def greedy_for_compact(source_pool):
    target_pool = []

    # calculating process
    calculating_round = 1
    while len(source_pool) > 1:
        print("Round " + str(calculating_round))

        results = []

        share_lock = multiprocessing.Manager().Lock()

        manager = multiprocessing.Manager()
        current_count = multiprocessing.Manager().Value("i", 0)
        total_count = multiprocessing.Manager().Value("i", int((len(source_pool) * (len(source_pool) - 1)) / 2))

        process_pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() - 1)

        greedy = []
        current_count = 0
        for row in range(len(source_pool) - 1):
            for col in range(row + 1, len(source_pool)):
                result = process_pool.apply_async(calculate, args=(source_pool, row, col, share_lock,
                                                                   current_count, total_count))
                results.append(result)

        process_pool.close()
        process_pool.join()

        calculator_map = {}
        overlap_matrix = [[0 for col in range(len(source_pool))] for row in range(len(source_pool) - 1)]

        for result in results:
            row, col, overlap = result.get()[0], result.get()[1], result.get()[2]
            fuser = result.get()[3]
            overlap_matrix[row][col] = overlap
            if overlap > 0:
                greedy.append([overlap, row, col])
                calculator_map[[row, col]] = fuser

        # no overlap
        for row in range(len(overlap_matrix)):
            if sum(overlap_matrix[row]) == 0:
                target_pool.append(source_pool[row])

        print()
        greedy.sort(reverse=True)

        available = [True for index in range(len(source_pool))]
        source_pool = []

        for fuse_information in greedy:
            index_1 = fuse_information[1]
            index_2 = fuse_information[2]
            if available[index_1] and available[index_2]:
                source_pool.append(calculator_map.get([index_1, index_2]).get_fuse_pattern_dna())
                available[index_1] = False
                available[index_2] = False

        print("Current source pool has " + str(len(source_pool)) + " Pattern DNA.")
        calculating_round += 1

    return target_pool


if __name__ == "__main__":
    # init process
    source = get_source_pool(read_file_path)
    # calculating process
    target = greedy_for_compact(source)
