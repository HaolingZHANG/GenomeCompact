"""
Name: inherent

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Function(s):
(1) Some inherent concepts.
"""

# double mapping of bases or degenerate bases.
base_mapping = [
    {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "M": "AC", "R": "AG", "W": "AT", "S": "CG", "Y": "CT", "K": "GT",
        "V": "ACG", "H": "ACT", "D": "AGT", "B": "CGT",
        "N": "ACGT"
    }, {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "AC": "M", "AG": "R", "AT": "W", "CG": "S", "CT": "Y", "GT": "K",
        "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B",
        "N": "ACGT"
    }
]

base_symbol = [
    "A", "C", "G", "T",
    "M", "R", "W", "S", "Y", "K",
    "V", "H", "D", "B",
    "N"
]

# mapping of complementary base
complementary = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "M": "K", "R": "Y", "W": "W", "S": "S", "Y": "R", "K": "M",
    "V": "B", "H": "D", "D": "H", "B": "V",
    "N": "N"
}

# mapping between normal proteins and codons.
protein_codon = {
    "A": ["GCN"],
    "C": ["TGY"],
    "D": ["GAY"],
    "E": ["GAR"],
    "F": ["TTY"],
    "G": ["GGN"],
    "H": ["CAY"],
    "I": ["ATH"],
    "K": ["AAR"],
    "L": ["TTR", "CTN"],
    "M": ["ATG"],
    "N": ["AAY"],
    "P": ["CCN"],
    "Q": ["CAR"],
    "R": ["CGN", "AGR"],
    "S": ["AGY", "TCN"],
    "T": ["ACN"],
    "V": ["GTN"],
    "W": ["TGG"],
    "Y": ["TAY"],
    "*": ["TAR", "TGA"]
}