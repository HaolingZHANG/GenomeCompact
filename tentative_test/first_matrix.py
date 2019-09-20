"""
Name: First Matrix

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Functions:
(1) Obtain the first overlap matrix to compare.
"""
import multiprocessing
from datetime import datetime

from methods.fuser import Fuser
from methods.fuser import init_indices
from methods.fuser import init_matched_group
from methods.dna import PatternDNA
from utils import accessor

read_path = "../data_set/Ecoli_K-1_MG1655.protein.fa"
write_path = "../output/first.csv"


def calculate_overlap(pool, row, col, lock, current, total, start):
    fuser = Fuser(init_indices(pool[row], pool[col]))
    overlap_length = fuser.calculate_overlap(init_matched_group(pool[row], pool[col]))

    lock.acquire()
    current.value += 1
    lock.release()
    print("\rdetect: " + str(current.value) + ", (" + str(total.value) +
          ") || (" + str(row) + ", " + str(col) + ") = " + str(overlap_length) +
          " || left time = " + str(datetime.now() - start), end=" ")

    return [row, col, overlap_length]


if __name__ == '__main__':
    start_time = datetime.now()

    source_pool = []
    indices, data_domains = accessor.read_information_by_fa(read_path)
    index = 0
    for index in range(len(data_domains)):
        source_pool.append(PatternDNA([index], protein=map(ord, data_domains[index])))
        index += 1

    results = []
    overlap_matrix = [[0 for col in range(len(source_pool))] for row in range(len(source_pool) - 1)]

    share_lock = multiprocessing.Manager().Lock()

    manager = multiprocessing.Manager()
    current_count = multiprocessing.Manager().Value("i", 0)
    total_count = multiprocessing.Manager().Value("i", int((len(source_pool) * (len(source_pool) - 1)) / 2))

    process_pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() - 1)

    for row in range(len(source_pool) - 1):
        for col in range(row + 1, len(source_pool)):
            result = process_pool.apply_async(calculate_overlap, args=(source_pool, row, col, share_lock,
                                                                       current_count, total_count, start_time))
            results.append(result)

    process_pool.close()
    process_pool.join()

    for result in results:
        overlap_matrix[result.get()[0]][result.get()[1]] = result.get()[2]

    print()
    end_time = datetime.now()
    print("The function run time is : %.03f seconds" % (end_time - start_time).seconds)

    # accessor.format_output("overlap matrix", overlap_matrix)
    with open(write_path, "w", encoding="utf-8") as save_file:
        for row in range(len(overlap_matrix)):
            row_data = str(overlap_matrix[row])[1: -1]
            save_file.write(row_data + "\n")
