"""
Name: inherent

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Function(s):
(1) Some inherent concepts.
"""
import numpy

# double mapping of bases or degenerate bases.
base_symbols = numpy.array([
    ["A", "C", "G", "T", "M",  "R",  "W",  "S",  "Y",  "K",  "V",   "H",   "D",   "B",   "N"],
    ["A", "C", "G", "T", "AC", "AG", "AT", "CG", "CT", "GT", "ACG", "ACT", "AGT", "CGT", "ACGT"]
])

# principle of complementary base pairing
pairing = numpy.array([
    ["A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N"],
    ["T", "G", "C", "A", "K", "Y", "W", "S", "R", "M", "B", "D", "H", "V", "N"],
])

# mapping between normal proteins and codons.
codons = numpy.array([
    ["A", "C", "D", "E", "F",
     "G", "H", "I", "K", "L",
     "L", "M", "N", "P", "Q",
     "R", "R", "S", "S",  "T",
     "V", "W", "Y", "*", "*"],
    ["GCN", "TGY", "GAY", "GAR", "TTY",
     "GGN", "CAY", "ATH", "AAR", "TTR",
     "CTN", "ATG", "AAY", "CCN", "CAR",
     "CGN", "AGR", "AGY", "TCN", "ACN",
     "GTN", "TGG", "TAY", "TAR", "TGA"]
])
