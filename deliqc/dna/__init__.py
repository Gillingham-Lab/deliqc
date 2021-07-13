import numpy as np
from .sanitation import is_valid, sanitise
from .codons import get_codon_coordinates

complements = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
}

bases = "ACGTN"
baseColour = {
    "G": np.array((31, 90, 173)) / 255,
    "C": np.array((80, 151, 250)) / 255,

    "A": np.array((173, 59, 45)) / 255,
    "T": np.array((250, 134, 120)) / 255,

    "N": (0.5, 0.5, 0.5),
}
baseToIndexDict = {bases[x]: x for x in range(len(bases))}


def get_nucleotide_colour_sequence(sequence):
    bar_colors = []

    for letter in sequence:
        try:
            bar_colors.append(baseColour[letter])
        except Exception:
            # print(letter)
            bar_colors.append(baseColour["N"])

    return bar_colors


def base2Index(x):
    return baseToIndexDict[x]


def reverse_complement(sequence: str) -> str:
    """
    Returns the reverse complement of a sequence
    :param sequence:
    :return:
    """

    return reverse(complement(sequence))


def reverse(sequence: str) -> str:
    """
    Returns the reverse of a sequence
    :param sequence:
    :return:
    """

    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Returns the complement of a sequence
    :param sequence:
    :return:
    """
    rc_seq = ""

    for base in sequence:
        rc_seq += complements[base]

    return rc_seq
