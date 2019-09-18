"""
Name: PatternDNA

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Functions:
(1) Get the max overlap of two Pattern DNAs;
(2) Fuse two Pattern DNAs.
"""

import sys
import copy

from methods.inherent import *
# from methods.dna import PatternDNA
# from methods.recoder import Recoder


# noinspection PyMethodMayBeStatic
class Fuser(object):

    def __init__(self, indices):
        """
        Initialize all variables

        ;:param indices: (fused) indices of two Pattern DNAs.
        """
        self._indices = indices
        self._major_strand = None
        self._minor_strand = None
        self._fuse_chooser = None
        self._start_position = -1
        self._overlap = 0

    # noinspection PyUnusedLocal
    def calculate_overlap(self, compare_matrix):
        """
        Get the maximum overlap of two Pattern DNAs.
        Ensure that overlapping parts can overlap completely, if not, overlap = 0.
        While calculating the maximum overlap, the smallest change (entropy) is also a considered factor (see replace).
        Currently, this method can not solve the protein with the intron.

        :param compare_matrix: one Pattern DNA.

        :return: the maximum overlap of the above two Pattern DNA.
        """
        min_replace = sys.maxsize

        for group in compare_matrix:
            major_strand, minor_strand = group[0], group[1]
            for start_position in range(len(major_strand[0])):
                replace = 0
                major_state = 0
                minor_state = 0
                major_fork_indices = numpy.where(major_strand[1] != 0)[0]
                minor_fork_indices = numpy.where(minor_strand[1] != 0)[0]

                if numpy.sum(major_fork_indices == start_position) > 0:
                    major_state = 3
                if numpy.sum(minor_fork_indices == 0) > 0:
                    minor_state = 3

                fork_chooser = numpy.array([[3 for index in range(int(len(major_fork_indices) / 3))],
                                            [3 for index in range(int(len(minor_fork_indices) / 3))]])

                for index in range(len(minor_strand[0])):
                    # print("start index = " + str(start_position + index) + "index = " + str(index))
                    # print("major_state = " + str(major_state) + ", minor_state = " + str(minor_state))
                    passes = self._get_current_passes(major_state, minor_state)

                    # print("init_pass_way = " + str(passes))

                    major_bases = self._get_current_bases(major_strand, start_position + index, major_state)
                    minor_bases = self._get_current_bases(minor_strand, index, minor_state)
                    same_flag = False

                    # print("major_bases = " + str(list(map(chr, major_bases))) +
                    #       ", minor_bases = " + str(list(map(chr, minor_bases))))

                    for major_index in range(len(major_bases)):
                        for minor_index in range(len(minor_bases)):
                            if passes[major_index][minor_index] \
                                    and is_same_base(major_bases[major_index], minor_bases[minor_index]):
                                same_flag = True
                            else:
                                passes[major_index][minor_index] = False
                                replace += 1

                    # print("result_pass_way = " + str(passes))

                    if not same_flag:
                        break

                    major_state, minor_state, fork_chooser = self._update_state(start_position, index, passes,
                                                                                major_state, minor_state,
                                                                                major_fork_indices, minor_fork_indices,
                                                                                fork_chooser)

                    if start_position + index == len(major_strand[0]) - 1 or index == len(minor_strand[0]) - 1:
                        if ((index + 1) > self._overlap) or ((index + 1) == self._overlap and min_replace > replace):
                            self._major_strand = copy.deepcopy(major_strand)
                            self._minor_strand = copy.deepcopy(minor_strand)
                            self._fuse_chooser = fork_chooser
                            self._start_position = start_position
                            self._overlap = index + 1
                            min_replace = replace
                        break

        return self._overlap

    def get_fuse_pattern_dna(self):
        """
        Calculate the fuse single DNA strand by two Pattern DNA and their recoder.

        :return: fuse single DNA strand.
        """
        if self._overlap > 0:
            # delete unusable path.
            choosers = self._fuse_chooser
            update_major_strand = self._update_strand(self._major_strand, choosers[0])
            update_minor_strand = self._update_strand(self._minor_strand, choosers[1])

            # fuse two DNA single strands.
            strand = self._fuse_strands(update_major_strand, update_minor_strand, self._start_position)

            return self._indices, strand

        return None

    def _update_strand(self, original_strand, fork_choosers):
        """
        Update the DNA single strand by its choosers

        :param original_strand:
        :param fork_choosers:

        :return: update Pattern DNA
        """
        strand = copy.deepcopy(original_strand)
        positions = numpy.where(original_strand[1] != 0)[0]
        index = 0
        for fork_chooser in fork_choosers:
            if fork_chooser < 3:
                if fork_chooser == 2:
                    strand[0][positions[index: index + 3]] = strand[1][positions[index: index + 3]]
                strand[1][positions[index: index + 3]] = 0
            index += 3

        return strand

    def _fuse_strands(self, major_strand, minor_strand, start_position):
        """
        Fuse two DNA single strands.

        :param major_strand:
        :param minor_strand:
        :param start_position:

        :return: fused DNA single strand.
        """

        # fuse prepositive strand
        strand = numpy.array([major_strand[0][:start_position], major_strand[1][:start_position]])

        # TODO consider the

        return strand

    def _get_current_bases(self, strand, position, state):
        """
        Get current available bases based on current position in the single strand.

        :param strand: the single strand in a Pattern DNA.
        :param position: comparison position of the single strand.
        :param state: state of current comparison position. 1 = only path 1, 2 = only path 2, 3 = both path.

        :return: current available bases.
        """

        current_bases = numpy.array([0, 0])

        if position <= len(strand[0]) - 1:
            if state == 0:
                current_bases[0] = strand[0][position]
            elif 1 <= state <= 2:
                current_bases[state - 1] = strand[0][position] if state <= 1 else strand[1][position]
            else:
                current_bases[0] = strand[0][position]
                current_bases[1] = strand[1][position]
        else:
            current_bases[0] = N

        return current_bases

    def _get_current_passes(self, major_state, minor_state):
        """
        Get the current passes based on the major state and minor state.

        :param major_state: major state of current comparison position. 1 = only path 1, 2 = only path 2, 3 = both path.
        :param minor_state: minor state of current comparison position. 1 = only path 1, 2 = only path 2, 3 = both path.

        :return: current passes.
        """
        passes = [[True, True], [True, True]]

        if 1 <= major_state <= 2:
            passes[(3 - major_state) - 1] = [False, False]
        if 1 <= minor_state <= 2:
            passes[0][(3 - minor_state) - 1] = False
            passes[1][(3 - minor_state) - 1] = False

        return passes

    def _update_state(self, start_position, index, passes, major_state, minor_state,
                        major_fork_indices, minor_fork_indices, fork_chooser):
        """
        Update the current major and minor state by the current position and current passes.

        :param start_position: start position in major strand.
        :param index: current reading position in minor strand.
        :param passes: current passes for two strand.
        :param major_state: original major state.
        :param minor_state:  original minor state.
        :param major_fork_indices: fork indices in major strand.
        :param minor_fork_indices: fork indices in minor strand.
        :param fork_chooser: original fork chooser.

        :return major_state:  update major state.
        :return minor_state: update minor state.
        :return fork_chooser: update fork chooser.
        """
        if major_state > 0:
            major_state = 2 * int(passes[1][0] or passes[1][1]) + int(passes[0][0] or passes[0][1])
        if minor_state > 0:
            minor_state = 2 * int(passes[0][1] or passes[1][1]) + int(passes[0][0] or passes[1][0])

        fork_situations = [
            numpy.where(major_fork_indices == start_position + index),
            numpy.where(major_fork_indices == start_position + index + 1),
            numpy.where(minor_fork_indices == index),
            numpy.where(minor_fork_indices == index + 1)
        ]

        if len(fork_situations[0]) > 0:
            fork_position = fork_situations[0][0]
            if fork_position % 3 == 2:
                fork_chooser[0][int(fork_position / 3)] = major_state
                major_state = 0
        if len(fork_situations[1]) > 0:
            fork_position = fork_situations[1][0]
            if fork_position % 3 == 0:
                major_state = 3

        if len(fork_situations[2]) > 0:
            fork_position = fork_situations[2][0]
            if fork_position % 3 == 2:
                fork_chooser[1][int(fork_position / 3)] = minor_state
                minor_state = 0
        if len(fork_situations[3]) > 0:
            fork_position = fork_situations[3][0]
            if fork_position % 3 == 0:
                minor_state = 3

        return major_state, minor_state, fork_chooser


def init_indices(dna_1, dna_2):
    """
    Obtain fuse indices.

    :param dna_1: one Pattern DNA.
    :param dna_2: another Pattern DNA.

    :return: the fuse indices.
    """
    return numpy.array(dna_1.get_indices() + dna_2.get_indices())


def init_matched_group(dna_1, dna_2):
    """
    Initialize DNA single strands to be compared.
    The major strand (as matrix[i][0]) describes the prepositive single DNA strand,
    and the minor strand (as matrix[i][1]) describes the postpositive single DNA strand.

    :param dna_1: one Pattern DNA.
    :param dna_2: another Pattern DNA.

    :return: compare matrix of the two above Pattern DNAs.
    """
    strand_1 = dna_1.get_strand()
    strand_2 = dna_2.get_strand(1)

    compare_matrix = []
    for index_1 in range(len(strand_1)):
        for index_2 in range(len(strand_2)):
            compare_matrix.append([copy.deepcopy(strand_1[index_1]), copy.deepcopy(strand_2[index_2])])
            compare_matrix.append([copy.deepcopy(strand_2[index_2]), copy.deepcopy(strand_1[index_1])])

    return compare_matrix


def is_same_base(base_1, base_2):
    """
    Equal judgement of two bases.
    Here, equality can mean that two characters of base are equal,
    or a degenerate base contains another base.

    :param base_1: one base used for comparison.
    :param base_2: another base used for comparison.

    :return: result of judgement.
    """

    # if base_1 is None or base_2 is None:
    if base_1 == 0 or base_2 == 0:
        return False

    # calculate same bases.
    if base_1 == base_2:
        return True

    # calculate inclusion.
    fused_values = get_intersection(base_1, base_2)

    # judge whether there is an intersection between the two base values.
    return len(fused_values) > 0


def fused_bases(base_1, base_2):
    """
    Get the intersecting base from two base.

    :param base_1: one base used for fusion.
    :param base_2: another base used for fusion.

    :return: base fused by two inumpyut base.
    """
    if base_1 == base_2:
        return base_1

    # calculate inclusion.
    fused_values = get_intersection(base_1, base_2)

    if len(fused_values) > 0:
        return simple_bases[0][numpy.where(simple_bases[1] == fused_values)][0]
    else:
        return 0


def get_intersection(base_1, base_2):
    """
    Get intersection of two base if possible.

    :param base_1: one base used for intersection.
    :param base_2: another base used for intersection.

    :return: intersection of two base.
    """

    base_1_values = detailed_bases[numpy.where(simple_bases == base_1)][0]
    base_2_values = detailed_bases[numpy.where(simple_bases == base_2)][0]
    intersection = numpy.intersect1d(base_1_values, base_2_values)

    return intersection
