from typing import List, Tuple


def get_codon_coordinates(dna) -> List[Tuple[int, int]]:
    """
    Returns a list of codon coordinates containing the start and the end of a codon.
    :param dna:
    :return:
    """
    coordinates = []

    in_codon = False
    start = None
    length = None

    for i in range(len(dna)):
        if in_codon is False:
            if dna[i] == "N":
                in_codon = True
                start = i
        else:
            if dna[i] != "N":
                coordinates.append((start, i))
                in_codon = False
                start = None

    return coordinates
