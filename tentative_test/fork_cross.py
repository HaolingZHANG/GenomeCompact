import numpy

from methods.dna import PatternDNA
from methods.fuser import Fuser
from methods.fuser import init_indices
from methods.fuser import init_matched_group

L = PatternDNA(["L"], dna_strand=numpy.array([["T", "T", "R"], ["C", "T", "N"]]))
R = PatternDNA(["R"], dna_strand=numpy.array([["C", "G", "N"], ["A", "G", "R"]]))
S = PatternDNA(["S"], dna_strand=numpy.array([["A", "G", "Y"], ["T", "C", "N"]]))
# stop codon
s = PatternDNA(["s"], dna_strand=numpy.array([["T", "A", "R"], ["T", "G", "A"]]))

fork_names = ["L", "R", "S", "s"]
forks = [L, R, S, s]

for index_1 in range(4):
    for index_2 in range(4):
        fuser = Fuser(init_indices(fork_names[index_2], forks[index_2]))
        overlap_length = fuser.calculate_overlap(init_matched_group(fork_names[index_2], forks[index_2]))
        print("maximum overlap = " + str(overlap_length))