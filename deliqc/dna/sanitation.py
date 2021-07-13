def is_valid(sequence: str) -> bool:
    """
    Checks if a given sequence contains only valid DNA.

    :param sequence:
    :return: True if sequence is valid.
    """
    elements = set("ACGTN")
    sequence_elements = set(sequence)

    if len(sequence) > 0 and len(elements | sequence_elements) == 5:
        return True
    else:
        return False


def sanitise(sequence):
    """
    Sanitises a given DNA sequence to be compatible with deliqc internals.

    :param sequence:
    :return: Sanitised sequence
    """
    return sequence.strip().upper()
