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

import numpy

from common.inherent import base_symbols
from common.dna import PatternDNA
from common.recoder import Recoder


# noinspection PyMethodMayBeStatic
class Fuser(object):

    def __init__(self):
        """
        Initialize all variables
        """
        self._indices = None
        self._major_strand = None
        self._minor_strand = None
        self._recoder = None
        self._overlap = 0

    def calculate_overlap(self, dna_1, dna_2):
        """
        Get the maximum overlap of two Pattern DNAs.
        Ensure that overlapping parts can overlap completely, if not, overlap = 0.
        While calculating the maximum overlap, the smallest change (entropy) is also a considered factor (see replace).
        Currently, this method can not solve the protein with the intron.

        :param dna_1: one Pattern DNA.
        :param dna_2: another Pattern DNA.

        :return: the maximum overlap of the above two Pattern DNA.
        """
        self._indices = dna_1.get_indices() + dna_2.get_indices()

        min_replace = sys.maxsize

        for group in self._init_matched_group(dna_1, dna_2):
            major_strand, minor_strand = group[0], group[1]
            for start_position in range(len(major_strand[0])):
                replace = 0
                recoder = Recoder(numpy.where(major_strand[1] != " ")[0],
                                  numpy.where(minor_strand[1] != " ")[0],
                                  start_position)

                for index in range(len(minor_strand[0])):
                    # print("index = " + str(index))
                    major_state, minor_state = recoder.get_current_state()
                    # print("major_state = " + str(major_state) + ", minor_state = " + str(minor_state))
                    passes = self._get_current_pass(major_state, minor_state)
                    # print("init_pass_way = " + str(passes))

                    major_bases = self._get_current_bases(major_strand, start_position + index, major_state)
                    minor_bases = self._get_current_bases(minor_strand, index, minor_state)
                    same_flag = False

                    # print("major_bases = " + str(major_bases) +
                    #       ", minor_bases = " + str(minor_bases))

                    for major_index in range(len(major_bases)):
                        for minor_index in range(len(minor_bases)):
                            if passes[major_index][minor_index] \
                                    and self.is_same_base(major_bases[major_index],
                                                          minor_bases[minor_index]):
                                same_flag = True
                            else:
                                passes[major_index][minor_index] = False
                                replace += 1

                    # print("result_pass_way = " + str(passes))

                    if not same_flag:
                        break

                    recoder.update(index, passes)
                    
                    if start_position + index == len(major_strand[0]) - 1 or index == len(minor_strand[0]) - 1:
                        if ((index + 1) > self._overlap) or ((index + 1) == self._overlap and min_replace > replace):
                            self._major_strand = copy.deepcopy(major_strand)
                            self._minor_strand = copy.deepcopy(minor_strand)
                            self._recoder = recoder
                            self._overlap = index + 1
                            min_replace = replace
                        break

    def _init_matched_group(self, dna_1, dna_2):
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

    def _get_current_bases(self, strand, position, state):
        """
        Get current available bases based on current position in the single strand.

        :param strand: the single strand in a Pattern DNA.
        :param position: comparison position of the single strand.
        :param state: state of current comparison position. 1 = only path 1, 2 = only path 2, 3 = both path.

        :return: current available bases.
        """
        if position <= len(strand[0]) - 1:
            current_bases = [" ", " "]
            if state == 0:
                current_bases[0] = strand[0][position]
            elif 1 <= state <= 2:
                current_bases[state - 1] = strand[0][position] if state <= 1 else strand[1][position]
            else:
                current_bases = [strand[0][position], strand[1][position]]

            return current_bases
        else:
            return ["N", " "]

    def _get_current_pass(self, major_state, minor_state):
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

    def get_fuse_pattern_dna(self):
        """
        Calculate the fuse single DNA strand by two Pattern DNA and their recoder.

        :return: fuse single DNA strand.
        """
        if self._overlap > 0:
            # delete unusable path.
            choosers = self._recoder.get_fuse_chooser()
            update_major_strand = self._update_strand(self._major_strand, choosers[0])
            update_minor_strand = self._update_strand(self._minor_strand, choosers[1])

            exit(111)
            # merge two DNA single strands.
            strand = self._merge_strands(update_major_strand, update_minor_strand, self._recoder.get_start_position())

            return PatternDNA(indices=self._indices, dna_strand=strand)

        return None

    def _update_strand(self, original_strand, fork_choosers):
        """
        Update the DNA single strand by its choosers

        :param original_strand:
        :param fork_choosers:

        :return: update Pattern DNA
        """
        # print("".join(original_strand[0]))
        # print("".join(original_strand[1]))
        strand = copy.deepcopy(original_strand)
        positions = numpy.where(original_strand[1] != " ")[0]
        index = 0
        for fork_chooser in fork_choosers:
            if fork_chooser < 3:
                if fork_chooser == 2:
                    strand[0][positions[index: index + 3]] = strand[1][positions[index: index + 3]]
                strand[1][positions[index: index + 3]] = " "
            index += 3

        # print("".join(strand[0]))
        # print("".join(strand[1]))
        return strand

    def _merge_strands(self, major_strand, minor_strand, start_position):
        """
        Merge two DNA single strands.

        :param major_strand:
        :param minor_strand:
        :param start_position:

        :return: merged
        """

        # fuse prepositive strand
        strand = [major_strand[0][:start_position], major_strand[1][:start_position]]
        # TODO consider the

        return strand

    def get_overlap_length(self):
        """
        Get the maximum overlap of the two above Pattern DNAs.

        :return: maximum overlap.
        """
        return self._overlap

    def is_same_base(self, base_1, base_2):
        """
        Equal judgement of two bases.
        Here, equality can mean that two characters of base are equal,
        or a degenerate base contains another base.

        :param base_1: one base used for comparison.
        :param base_2: another base used for comparison.

        :return: the result of judgement.
        """
        print(str(base_1) + ", " + str(base_1))

        if base_1 == " " or base_2 == " ":
            return False

        # calculate same bases.
        if base_1 == base_2:
            return True

        # calculate inclusion.
        fused_values = self._get_intersection(base_1, base_2)

        # judge whether there is an intersection between the two base values.
        return len(fused_values) > 0

    def fused_bases(self, base_1, base_2):
        """
        Get the intersecting base from two base.

        :param base_1: one base used for fusion.
        :param base_2: another base used for fusion.

        :return: the base fused by two input base.
        """
        if base_1 == base_2:
            return base_1

        # calculate inclusion.
        fused_values = self._get_intersection(base_1, base_2)

        if len(fused_values) > 0:
            return base_symbols[0][numpy.where(base_symbols[1] == fused_values)][0]
        else:
            return " "

    def _get_intersection(self, base_1, base_2):
        base_1_values = numpy.array(list(base_symbols[1][numpy.where(base_symbols[0] == base_1)][0]))
        base_2_values = numpy.array(list(base_symbols[1][numpy.where(base_symbols[0] == base_2)][0]))
        intersection = numpy.intersect1d(base_1_values, base_2_values)
        return intersection
