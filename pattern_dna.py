"""
Name: PatternDNA

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Functions:
(1) Initiate Pattern DNA from a DNA single strand or a protein;
(2) Get Information of the created Pattern DNA.
"""

from inherent import *


# noinspection PyMethodMayBeStatic
class PatternDNA(object):

    def __init__(self, indices, protein=None, dna_strand=None):
        """
        Initiate the pattern DNA by protein or DNA single strand.
        A binary choice, input protein or dna single strand.

        :param protein: the input protein.
        :param dna_strand: the obtain DNA single strand.
        """
        if protein is None and dna_strand is None:
            raise ValueError("Pattern DNA " + str(indices) + " is Error, because we don't obtain value.")
        elif protein is not None and dna_strand is not None:
            raise ValueError("Pattern DNA " + str(indices) + " is Error, because we don't which value can be selected.")

        self._indices = indices

        if protein is not None:
            self._template_strand = self._get_template_strand(protein)
            self._complementary_strand = self._get_complementary_strand(self._template_strand)
        else:
            self._template_strand = dna_strand
            self._complementary_strand = self._get_complementary_strand(dna_strand)

    def _get_template_strand(self, protein):
        """
        Get the corresponding template DNA strand from the protein.
        The data structure of DNA single strand describes as {
            "s": strand,
            "m_i": multi_indices,
            "m_b": multi_bases,
            "l": len(strand)}.

        :param protein: the obtained protein.

        :return: corresponding information of the template DNA strand.
        """
        strand = []
        multi_indices = []
        multi_bases = {}

        for index in range(len(protein)):
            codon_chooser = protein_codon.get(protein[index])
            for nucleotide in codon_chooser[0]:
                strand.append(nucleotide)
            if len(codon_chooser) > 1:
                for codon_index in range(3):
                    multi_index = len(strand) - 3 + codon_index
                    multi_indices.append(multi_index)
                    multi_bases[multi_index] = codon_chooser[1][codon_index]

        information = {
            "s": strand,
            "m_i": multi_indices,
            "m_b": multi_bases,
            "l": len(strand)
        }

        return information

    def _get_complementary_strand(self, input_strand):
        """
        Get the corresponding complementary DNA strand from the normal DNA strand.
        The data structure of DNA single strand describes as {
            "s": strand,
            "m_i": multi_indices,
            "m_b": multi_bases,
            "l": len(strand)}.

        :param input_strand: normal strand of the pattern DNA.

        :return: corresponding information of complementary DNA strand.
        """
        strand = []
        multi_indices = []
        multi_bases = {}

        for base in input_strand.get("s"):
            strand.append(complementary.get(base))
        strand = strand[::-1]

        for key, multi_base in input_strand.get("m_b").items():
            complementary_key = len(strand) - 1 - int(key)
            multi_indices.append(complementary_key)
            multi_bases[complementary_key] = complementary.get(multi_base)

        information = {
            "s": strand,
            "m_i": sorted(multi_indices),
            "m_b": multi_bases,
            "l": len(strand)
        }

        return information

    def get_template_strand(self):
        """
        Get the template strand.

        :return: template strand.
        """
        return self._template_strand

    def get_complementary_strand(self):
        """
        Get the complementary strand.

        :return: complementary strand.
        """
        return self._complementary_strand

    def get_dna_information(self, start_position):
        """
        Get the intuitive images of Pattern DNA.

        :param start_position: start position of reading the Pattern DNA.

        :return: intuitive strand information of this Pattern DNA.
        """
        information = []
        template_multi = ""
        complementary_multi = ""

        for index in range(self._template_strand.get("l")):
            if index in self._template_strand.get("m_i"):
                template_multi += self._template_strand.get("m_b").get(index)
            else:
                template_multi += " "
        information.append(["".join(self._template_strand.get("s")[start_position:]),
                            template_multi[start_position:]])

        for index in range(self._complementary_strand.get("l")):
            if index in self._complementary_strand.get("m_i"):
                complementary_multi += self._complementary_strand.get("m_b").get(index)
            else:
                complementary_multi += " "
        information.append(["".join(self._complementary_strand.get("s")[start_position:]),
                            complementary_multi[start_position:]])

        return information

    def get_indices(self):
        """
        Get the protein indices of this Pattern DNA.

        :return: protein indices of this Pattern DNA.
        """
        return self._indices

    def __str__(self):
        """
        Get the whole intuitive information of this Pattern DNA.

        :return: whole intuitive information of this Pattern DNA.
        """
        information = self.get_dna_information(0)

        return "protein indices = " + str(self._indices) + ": \n" + \
               "t~ strand = " + str(information[0][0]) + "\n" + \
               "            " + str(information[0][1]) + "\n" + \
               "c~ strand = " + str(information[1][0]) + "\n" + \
               "            " + str(information[1][1]) + "\n"
