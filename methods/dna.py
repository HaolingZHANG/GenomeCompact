"""
Name: PatternDNA

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Functions:
(1) Initiate Pattern DNA from a DNA single strand or a protein;
(2) Get Information of the created Pattern DNA.
"""
import copy

from methods.inherent import *


# noinspection PyMethodMayBeStatic
class PatternDNA(object):

    def __init__(self, indices, protein=None, dna_strand=None):
        """
        Initiate the pattern DNA by protein or DNA single strand.
        A binary choice, inumpyut protein or dna single strand.

        :param protein: the inumpyut protein (int[:]).
        :param dna_strand: the obtain DNA single strand (int[:]).
        """
        if protein is None and dna_strand is None:
            raise ValueError("Pattern DNA " + str(indices) + " is Error, because the value cannot be obtained.")
        elif protein is not None and dna_strand is not None:
            raise ValueError("Pattern DNA " + str(indices) + " is Error, because the value cannot be selected.")

        self._indices = indices

        if protein is not None:
            self._t_strand = self._obtain_t_strand(protein)
            self._c_strand = self._obtain_c_strand(self._t_strand)
        else:
            self._t_strand = copy.deepcopy(dna_strand)
            self._c_strand = self._obtain_c_strand(dna_strand)

    def _obtain_t_strand(self, protein):
        """
        Get the corresponding template DNA strand from the protein.

        :param protein: the obtained protein.

        :return: corresponding information of the template DNA strand.
        """

        strand = [[], []]
        for amino_acid in protein:
            codon_chooser = base_maps[numpy.where(codons == amino_acid)]
            for index in range(3):
                strand[0].append(codon_chooser[0][index])
                if len(codon_chooser) > 1:
                    strand[1].append(codon_chooser[1][index])
                else:
                    strand[1].append(0)

        return numpy.array(strand)

    def _obtain_c_strand(self, template_strand):
        """
        Get the corresponding complementary DNA strand from the normal DNA strand.

        :param template_strand: normal strand of the pattern DNA.

        :return: corresponding information of complementary DNA strand.
        """
        input_strand = numpy.transpose(copy.deepcopy(template_strand))
        strand = [[], []]
        for bases in input_strand:
            strand[0].insert(0, c_pairing[numpy.where(t_pairing == bases[0])][0])
            if bases[1] != 0:
                strand[1].insert(0, c_pairing[numpy.where(t_pairing == bases[1])][0])
            else:
                strand[1].insert(0, 0)

        return numpy.array(strand)

    def get_strand(self, strand_type=0):
        """
        Get the strand by strand type.
        1 is the template strand; -1 is the complementary strand; the others is the whole strand.

        :return: strand.
        """
        if strand_type == 1:
            return [self._t_strand]
        elif strand_type == -1:
            return [self._c_strand]
        else:
            return [self._t_strand, self._c_strand]

    def get_dna_fragment(self, start_position=0, stop_position=0):
        """
        Get the intuitive images of Pattern DNA.

        :param start_position: start position of reading the Pattern DNA.
        :param stop_position: stop position of reading the Pattern DNA. (positive number)

        :return: intuitive strand information of this Pattern DNA.
        """

        if start_position != 0:
            information = [["".join([chr(data) for data in self._t_strand[0]][start_position: ]),
                            "".join([chr(data) for data in self._t_strand[1]][start_position: ])],
                           ["".join([chr(data) for data in self._c_strand[0]][start_position: ]),
                            "".join([chr(data) for data in self._c_strand[1]][start_position: ])]]
        else:
            information = [["".join([chr(data) for data in self._t_strand[0]][start_position:]),
                            "".join([chr(data) for data in self._t_strand[1]][start_position:])],
                           ["".join([chr(data) for data in self._c_strand[0]][start_position:]),
                            "".join([chr(data) for data in self._c_strand[1]][start_position:])]]

        return information

    def get_indices(self):
        """
        Get the protein indices of this Pattern DNA.

        :return: protein indices of this Pattern DNA.
        """
        return self._indices

    def get_all_proteins(self):
        pass

    def __str__(self):
        """
        Get the whole intuitive information of this Pattern DNA.

        :return: whole intuitive information of this Pattern DNA.
        """
        information = self.get_dna_fragment()

        return "protein indices = " + str(self._indices) + ": \n" + \
               "t~ strand = " + str(information[0][0]) + "\n" + \
               "            " + str(information[0][1]) + "\n" + \
               "c~ strand = " + str(information[1][0]) + "\n" + \
               "            " + str(information[1][1]) + "\n"
