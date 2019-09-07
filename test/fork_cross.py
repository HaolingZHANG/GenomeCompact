import numpy

from common.dna import PatternDNA
from common.fuser import Fuser

L = PatternDNA(["L"], dna_strand=numpy.array([["T", "T", "R"], ["C", "T", "N"]]))
R = PatternDNA(["R"], dna_strand=numpy.array([["C", "G", "N"], ["A", "G", "R"]]))
S = PatternDNA(["S"], dna_strand=numpy.array([["A", "G", "Y"], ["T", "C", "N"]]))
# stop codon
s = PatternDNA(["s"], dna_strand=numpy.array([["T", "A", "R"], ["T", "G", "A"]]))

fork_names = ["L", "R", "S", "s"]
forks = [L, R, S, s]

for index_1 in range(4):
    for index_2 in range(4):
        fuser = Fuser()
        print(fork_names[index_1] + ", " + fork_names[index_2])
        fuser.calculate_overlap(forks[index_1], forks[index_2])
        print("maximum overlap = " + str(fuser.get_overlap_length()))
        print()