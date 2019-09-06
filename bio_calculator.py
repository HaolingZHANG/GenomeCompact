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

from numba import jit
from pattern_dna import PatternDNA

from bio_operator import *


# noinspection PyMethodMayBeStatic,PyUnusedLocal
class Calculator():

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
            overlap, best_start_position = 0, -1

            for start_position in range(major_strand.get("l")):
                replace = 0

                recoder = Calculator.Recoder(major_strand.get("m_i"), minor_strand.get("m_i"), start_position)

                for index in range(minor_strand.get("l")):
                    if start_position + index < major_strand.get("l") - 1 and index < minor_strand.get("l") - 1:
                        # print("index = " + str(index))
                        major_state, minor_state = recoder.get_current_state()
                        # print("major_state = " + str(major_state) + ", minor_state = " + str(minor_state))
                        passes = self._get_current_pass(major_state, minor_state)
                        # print("init_pass_way = " + str(passes))

                        current_major_bases = self._get_current_bases(major_strand, start_position + index, major_state)
                        current_minor_bases = self._get_current_bases(minor_strand, index, minor_state)
                        same_flag = False

                        # print("current_major_bases = " + str(current_major_bases) +
                        #       ", current_minor_bases = " + str(current_minor_bases))

                        for major_index in range(len(current_major_bases)):
                            for minor_index in range(len(current_minor_bases)):
                                if passes[major_index][minor_index] and is_same_base(current_major_bases[major_index],
                                                                                     current_minor_bases[minor_index]):
                                    same_flag = True
                                else:
                                    passes[major_index][minor_index] = False
                                    replace += 1

                        # print("result_pass_way = " + str(passes))

                        if not same_flag:
                            break

                        recoder.update(index, passes)
                    else:
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
        strand_1 = [dna_1.get_template_strand(), dna_1.get_complementary_strand()]
        strand_2 = [dna_2.get_template_strand()]

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
        current_bases = [None, None]
        if state == 0:
            current_bases[0] = strand.get("s")[position]
        elif 1 <= state <= 2:
            current_bases[state - 1] = strand.get("s")[position] if state <= 1 else strand.get("m_b").get(position)
        else:
            current_bases[0] = strand.get("s")[position]
            current_bases[1] = strand.get("m_b").get(position)

        return current_bases

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

    def get_fuse_dna(self):
        """
        Calculate the fuse single DNA strand by two Pattern DNA and their recoder.

        :return: fuse single DNA strand.
        """
        if self._overlap > 0:
            # delete unusable path.
            choosers = self._recoder.get_fuse_chooser()
            major_strand = self._update_strand(self._major_strand, choosers[0])
            minor_strand = self._update_strand(self._minor_strand, choosers[1])

            # merge two DNA single strands.
            strand = self._merge_strands(major_strand, minor_strand, self._recoder.get_start_position())

            return PatternDNA(indices=self._indices, dna_strand=strand)

        return None

    def _update_strand(self, input_strand, choosers):
        """
        Update the DNA single strand by its choosers

        :param input_strand:
        :param choosers:

        :return: update Pattern DNA
        """
        strand = copy.deepcopy(input_strand.get("s"))
        multi_indices = copy.deepcopy(input_strand.get("m_i"))
        multi_bases = copy.deepcopy(input_strand.get("m_b"))

        fork_index = 0
        for chooser in choosers:
            if chooser < 3:
                if chooser == 2:
                    for index in range(fork_index, fork_index + 3):
                        strand[multi_indices[index]] = multi_bases.get(multi_indices[index])

                delete_indices = []
                for index in range(fork_index, fork_index + 3):
                    delete_indices.append(multi_indices[index])
                    del multi_bases[multi_indices[index]]
                for delete_index in delete_indices:
                    multi_indices.remove(delete_index)
            fork_index += 3

        information = {
            "s": strand,
            "m_i": sorted(multi_indices),
            "m_b": multi_bases,
            "l": len(strand)
        }

        return information

    def _merge_strands(self, major_strand, minor_strand, start_position):
        """
        Merge two DNA single strands.

        :param major_strand:
        :param minor_strand:
        :param start_position:

        :return: merged
        """
        strand = []
        multi_indices = []
        multi_bases = {}

        # fuse prepositive strand
        for major_index in range(start_position):
            strand.append(major_strand.get("s")[major_index])
            if major_index in major_strand.get("m_i"):
                multi_indices.append(major_index)
                multi_bases[major_index] = major_strand.get("m_b").get(major_index)

        # fuse overlap strand
        for fuse_index in range(self._overlap):
            current_major_bases = major_strand.get("s")[fuse_index + start_position]
            if (fuse_index + start_position) in major_strand.get("m_i"):
                current_major_bases.append(major_strand.get("m_b").get(fuse_index + start_position))
            current_minor_bases = minor_strand.get("s")[fuse_index]
            if fuse_index in minor_strand.get("m_i"):
                current_minor_bases.append(minor_strand.get("m_i").get(fuse_index))
            # TODO
            if len(current_major_bases) * len(current_minor_bases) == 1:
                strand.append(fused_bases(current_major_bases[0], current_minor_bases[0]))
            elif len(current_major_bases) * len(current_minor_bases) == 2:
                current_bases = []
                for major_base in current_major_bases:
                    for minor_base in current_minor_bases:
                        current_bases.append(fused_bases(major_base, minor_base))
                strand.append(current_bases[0])
                multi_indices.append(len(strand) - 1)
                multi_bases[len(strand) - 1] = current_bases[1]
            else:
                pass

        # fuse postpositive strand
        for minor_index in range(minor_strand.get("l") - self._overlap):
            strand.append(minor_strand.get("s")[minor_index + self._overlap])
            if (minor_index + self._overlap) in minor_strand.get("m_i"):
                multi_indices.append(minor_index + major_strand.get("l"))
                multi_bases[minor_index + major_strand.get("l")] = major_strand.get("m_b").get(minor_index +
                                                                                               major_strand.get("l"))

        information = {
            "s": strand,
            "m_i": sorted(multi_indices),
            "m_b": multi_bases,
            "l": len(strand)
        }

        return information

    def get_overlap_length(self):
        """
        Get the maximum overlap of the two above Pattern DNAs.

        :return: maximum overlap.
        """
        return self._overlap

    # noinspection PyUnusedLocal
    class Recoder(object):

        def __init__(self, major_fork_indices, minor_fork_indices, start_position):
            """
            Initialize recoder for save process variables in the process of calculating maximum overlap.

            :param major_fork_indices: index list with fork of major strand.
            :param minor_fork_indices: index list with fork of minor strand.
            :param start_position: the start position of major strand.
            """
            self._major_fork_indices = copy.deepcopy(major_fork_indices)
            self._minor_fork_indices = copy.deepcopy(minor_fork_indices)
            self._major_state = 0
            self._minor_state = 0

            self._start_position = start_position

            if start_position in major_fork_indices:
                self._major_state = 3

            if 0 in minor_fork_indices:
                self._minor_state = 3

            self._major_fork_chooser = [3 for index in range(int(len(major_fork_indices) / 3))]
            self._minor_fork_chooser = [3 for index in range(int(len(minor_fork_indices) / 3))]

        def update(self, current_position, passes):
            """
            Update the current major and minor state by the current position and current passes.

            :param current_position: current reading position in minor strand/
            :param passes: current passes for two strand.
            """
            if self._major_state > 0:
                self._major_state = 0
                if passes[0][0] or passes[0][1]:
                    self._major_state += 1
                if passes[1][0] or passes[1][1]:
                    self._major_state += 2
            if self._minor_state > 0:
                self._minor_state = 0
                if passes[0][0] or passes[1][0]:
                    self._minor_state += 1
                if passes[0][1] or passes[1][1]:
                    self._minor_state += 2

            if (self._start_position + current_position) in self._major_fork_indices:
                fork_position = self._major_fork_indices.index(self._start_position + current_position)
                if fork_position % 3 == 2:
                    self._major_fork_chooser[int(fork_position / 3)] = self._major_state
                    self._major_state = 0
            if (self._start_position + current_position + 1) in self._major_fork_indices:
                fork_position = self._major_fork_indices.index(self._start_position + current_position + 1)
                if fork_position % 3 == 0:
                    self._major_state = 3
            if current_position == 0 and self._start_position in self._major_fork_indices:
                self._major_state = 3

            if current_position in self._minor_fork_indices:
                fork_position = self._minor_fork_indices.index(current_position)
                if fork_position % 3 == 2:
                    self._minor_fork_indices[int(fork_position / 3)] = self._minor_state
                    self._minor_state = 0
            if (current_position + 1) in self._minor_fork_indices:
                fork_position = self._minor_fork_indices.index(current_position + 1)
                if fork_position % 3 == 0:
                    self._minor_state = 3

        def get_current_state(self):
            """
            Get current major state and minor state.

            :return: current major state and minor state.
            """
            return self._major_state, self._minor_state

        def get_start_position(self):
            """
            Get the start position of this recoder.

            :return: start position.
            """
            return self._start_position

        def get_fuse_chooser(self):
            """
            Get the fuse chooser.

            :return: fuse chooser of major fork and minor fork.
            """
            return [self._major_fork_chooser, self._minor_fork_chooser]
