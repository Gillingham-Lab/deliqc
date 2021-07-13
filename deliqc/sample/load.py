from collections import defaultdict
import numpy as np
from typing import List
from multiprocessing import Pool

from .worker import worker, WorkerData
from .reader import read_lines
from .result import AlignedReadResult, FailedReadResult


def file_reader(fn1, fn2, num_reads, metadata):
    read = 0
    for i, (r1, r2) in enumerate(zip(read_lines(fn1, num_reads), read_lines(fn2, num_reads))):
        if read >= num_reads:
            break

        yield i, (r1, r2), metadata

        read += 1


def load_sample(
        reference: str,
        fn1: str,
        fn2: str,
        max_reads: int = 2000,
        max_pair_mismatches: int = 0,
        max_point_mutations: int = 5,
        split_on_codon: int = 0,
        codons: List[str] = [],
        threads: int = 8,
        batches: int = 100,
):
    read_length = len(reference)
    # Array to keep track of mismatch count
    total_mismatches = np.zeros((read_length, 1 + len(codons)))
    # Array to keep track of insertions (0) and deletions (1)
    total_indels = np.zeros((read_length, 2, 1 + len(codons)))
    # Array to keep track of the target of a point mutation
    mutation_targets = np.zeros((read_length, 5, 1 + len(codons)))
    # Array to keep track of the inserted base
    insertion_targets = np.zeros((read_length, 5, 1 + len(codons)))

    codon_counts = defaultdict(int)
    count_mismatched_pairs = 0
    count_too_many_point_mutations = 0
    count_aligned_reads = 0

    worker_data = WorkerData(
        reference=reference,
        max_pair_mismatches=max_pair_mismatches,
        prefer_mate=False,
        max_point_mutations=max_point_mutations,
        split_on_codon=split_on_codon,
        codons=codons,
    )

    reads_worked = 0
    with Pool(threads) as pool:
        # Prepare the file reader to seed workers
        r = file_reader(fn1, fn2, max_reads, worker_data)

        for result in pool.imap_unordered(worker, r, batches):
            reads_worked += 1

            # Check what kind of result we got

            if isinstance(result, AlignedReadResult):
                # Successfully aligned result
                count_aligned_reads += 1

                # Get codon index
                if result.codon in codons:
                    codon_index = codons.index(result.codon) + 1
                    codon_counts[result.codon] += 1
                else:
                    codon_index = 0

                # Count up
                total_mismatches[:, 0] += result.mismatches
                mutation_targets[:, :, 0] += result.mutationTarget
                total_indels[:, :, 0] += result.indels
                insertion_targets[:, :, 0] += result.insertionTarget

                # Count codon specific result
                if codon_index > 0:
                    total_mismatches[:, codon_index] += result.mismatches
                    mutation_targets[:, :, codon_index] += result.mutationTarget
                    total_indels[:, :, codon_index] += result.indels
                    insertion_targets[:, :, codon_index] += result.insertionTarget
            elif isinstance(result, FailedReadResult):
                # Failed result
                if result.reason == FailedReadResult.tooManyMismatchesBetweenPair:
                    count_mismatched_pairs += 1
                elif result.reason == FailedReadResult.tooManyMismatchesToReference:
                    count_too_many_point_mutations += 1

    # Normalise all
    total_mismatches[:, 0] /= count_aligned_reads
    mutation_targets[:, :, 0] /= count_aligned_reads
    total_indels[:, :, 0] /= count_aligned_reads
    insertion_targets[:, :, 0] /= count_aligned_reads

    # And normalise for each codon
    for i in range(len(codons)):
        codon = codons[i]

        total_mismatches[:, i+1] /= codon_counts[codon]
        mutation_targets[:, :, i+1] /= codon_counts[codon]
        total_indels[:, :, i+1] /= codon_counts[codon]
        insertion_targets[:, :, i+1] /= codon_counts[codon]

    return {
        "reads": reads_worked,
        "mismatchedPairs": count_mismatched_pairs,
        "tooManyPointMutations": count_too_many_point_mutations,
        "mismatches": total_mismatches,
        "indels": total_indels,
        "alignedReads": count_aligned_reads,
        "alignedReadsPerCodon": codon_counts,
        "mutationTarget": mutation_targets,
        "insertionTarget": insertion_targets,
    }

