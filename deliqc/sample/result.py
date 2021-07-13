import numpy as np


class BaseResult:
    pass


class AlignedReadResult(BaseResult):
    def __init__(self, L: int):
        self.mismatches = np.zeros(L)
        self.mutationTarget = np.zeros((L, 5))
        self.insertionTarget = np.zeros((L, 5))
        self.indels = np.zeros((L, 2))

        self.pairMismatches = 0
        self.codon = None
        self.codons = []

    @property
    def mismatchCount(self):
        return self.mismatches.sum()


class FailedReadResult(BaseResult):
    tooManyMismatchesBetweenPair = 0
    tooManyMismatchesToReference = 1

    def __init__(self, reason):
        self.reason = reason
