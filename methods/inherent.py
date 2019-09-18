"""
Name: inherent

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Function(s):
(1) Some inherent concepts.
"""
import numpy

# mapping of integer and char
A = 65  # ord('A')
B = 66  # ord('B')
C = 67  # ord('C')
D = 68  # ord('D')
E = 69  # ord('E')
F = 70  # ord('F')
G = 71  # ord('G')
H = 72  # ord('H')
I = 73  # ord('I')
K = 75  # ord('K')
L = 76  # ord('L')
M = 77  # ord('M')
N = 78  # ord('N')
P = 80  # ord('P')
Q = 81  # ord('Q')
R = 82  # ord('R')
S = 83  # ord('S')
T = 84  # ord('T')
V = 86  # ord('V')
W = 87  # ord('W')
Y = 89  # ord('Y')
e = 42  # ord('*')

# double mapping of bases or degenerate bases.
simple_bases = numpy.array([
    A, C, G, T,
    M, R, W, S, Y, K,
    V, H, D, B,
    N
])

detailed_bases = numpy.array([
    [A], [C], [G], [T],
    [A, C], [A, G], [A, T], [C, G], [C, T], [G, T],
    [A, C, G], [A, C, T], [A, G, T], [C, G, T],
    [A, C, G, T]
])


# principle of complementary base pairing
t_pairing = numpy.array([
    A, C, G, T, M, R, W, S, Y, K, V, H, D, B, N
])
c_pairing = numpy.array([
    T, G, C, A, K, Y, W, S, R, M, B, D, H, V, N
])

# mapping between normal proteins and codons.
codons = numpy.array([
    A, C, D, E, F,
    G, H, I, K, L,
    L, M, N, P, Q,
    R, R, S, S, T,
    V, W, Y, e, e
])
base_maps = numpy.array([
    [G, C, N], [T, G, Y], [G, A, Y], [G, A, R], [T, T, Y],
    [G, G, N], [C, A, Y], [A, T, H], [A, A, R], [T, T, R],
    [C, T, N], [A, T, G], [A, A, Y], [C, C, N], [C, A, R],
    [C, G, N], [A, G, R], [A, G, Y], [T, C, N], [A, C, N],
    [G, T, N], [T, G, G], [T, A, Y], [T, A, R], [T, G, A]
])
