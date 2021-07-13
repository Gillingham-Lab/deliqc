from typing import Tuple
from deliqc import dna
from deliqc.sample.alignment import _simple_analyzer, _advanced_analyzer
from deliqc.sample.result import FailedReadResult, AlignedReadResult
from deliqc.sample.errors import TooManyMismatchesBetweenPair, TooManyMismatchesToReference, WeirdReads


class WorkerData:
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.kwargs["coordinates"] = dna.get_codon_coordinates(kwargs["reference"])
        self.kwargs["read_length"] = len(kwargs["reference"])


def worker(data: Tuple[int, Tuple[str, str], WorkerData]):
    read_number, (r1, r2), metadata = data
    kwargs = metadata.kwargs

    simple_failed = False
    advanced_failed = False
    too_many_mismatches_between_pair = False
    too_many_mismatches_to_reference = False

    ref_mismatches = 0
    pair_mismatches = 0
    read_indel_target = None

    try:
        result = _simple_analyzer(r1, r2, **kwargs)
    except TooManyMismatchesBetweenPair:
        simple_failed = True
        too_many_mismatches_between_pair = True
    except TooManyMismatchesToReference:
        simple_failed = True
        too_many_mismatches_to_reference = True

    if simple_failed is True:
        try:
            result = _advanced_analyzer(r1, r2, **kwargs)
        except TooManyMismatchesBetweenPair:
            advanced_failed = True
            too_many_mismatches_between_pair = True
        except TooManyMismatchesToReference:
            advanced_failed = True
            too_many_mismatches_to_reference = True
        except WeirdReads:
            advanced_failed = True
        else:
            too_many_mismatches_between_pair = False
            too_many_mismatches_to_reference = False

    if advanced_failed is True:
        if too_many_mismatches_between_pair:
            result = FailedReadResult(FailedReadResult.tooManyMismatchesBetweenPair)

        if too_many_mismatches_to_reference:
            result = FailedReadResult(FailedReadResult.tooManyMismatchesToReference)

    return result
