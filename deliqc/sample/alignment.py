import numpy as np

from deliqc import dna
from deliqc.sample.result import AlignedReadResult, FailedReadResult
from deliqc.sample.errors import *


def _simple_analyzer(
        r1,
        r2,
        read_length,
        reference,
        max_pair_mismatches,
        max_point_mutations,
        prefer_mate,
        split_on_codon=False,
        codons=[],
        coordinates=None,
):
    """ This simple analyzer only looks at point mutations. """
    r1 = dna.reverse_complement(r1[:read_length])
    r2 = r2[:read_length]

    codon = None

    # Count mismatches between the two mates
    pairMismatches = 0
    realRead = ""
    for nt1, nt2 in zip(r1, r2):
        if nt1 != nt2:
            pairMismatches += 1

            realRead += nt1 if prefer_mate is False else nt2
        else:
            realRead += nt1

    if pairMismatches > max_pair_mismatches:
        raise TooManyMismatchesBetweenPair

    # Count the amount of mismatches compared to template
    result = AlignedReadResult(read_length)
    result.pairMismatches = pairMismatches

    # readMismatches = np.zeros(readLength)
    # readMutationTarget = np.zeros((readLength, 5))

    for i in range(read_length):
        try:
            nt_ref = reference[i]
            nt_read = realRead[i]
        except IndexError:
            print("Read too short?", realRead, reference)
            break

        # Continue if matching
        if nt_ref == nt_read:
            continue

        # Continue if reference contains N (=codon)
        if nt_ref == "N":
            continue

        result.mismatches[i] += 1
        result.mutationTarget[i, dna.base2Index(nt_read)] += 1

    if result.mismatchCount > max_point_mutations:
        raise TooManyMismatchesToReference

    # Split on codon if given
    codon = None
    if split_on_codon > 0:
        if split_on_codon < len(coordinates):
            try:
                x, y = coordinates[split_on_codon]
                codonRead = realRead[x:y]

                if codonRead in codons:
                    codon = codonRead
            except TypeError:
                print(type(x), x, type(y), y)
        else:
            print("Codon number is beyond limits.")

    result.codon = codon

    # Save other codons, too
    codons = []
    for x, y in coordinates:
        codons.append(realRead[x:y])
    result.codons = codon

    return result


def _match_score(x, y):
    if x == "N":
        return 1
    elif y == x:
        return 5
    else:
        return -2


def _advanced_analyzer(
        r1,
        r2,
        read_length,
        reference,
        max_pair_mismatches,
        max_point_mutations,
        prefer_mate,
        split_on_codon=False,
        codons=[],
        coordinates=None,
):
    """ This advanced analyzer aligns both read to the reference. Thus, it looks at point mutations AND indels."""
    from Bio import pairwise2

    # Invert read 1 (reverse read)
    r1 = dna.reverse_complement(r1)

    # Canonicalize both reads
    pair_a = pairwise2.align.localcs(r1, r2, _match_score, -10, -1, one_alignment_only=True)[0]

    ar1 = pair_a[0]
    ar2 = pair_a[1]

    # We first clip away overhangs coming from adapter sequences
    # Check if only ar1 starts with "-" and ar2 ends with "-"
    if ar1[0] == "-" and ar1[-1] != "-" and ar2[0] != "-" and ar2[-1] == "-":
        left = 0
        for x in ar1:
            if x == "-":
                left += 1
            else:
                break

        right = 0
        for x in ar2[::-1]:
            if x == "-":
                right += 1
            else:
                break

        length = len(ar1) - left - right

        ar1 = ar1[left:left + length]
        ar2 = ar2[left:left + length]
    # Check the inverse.
    elif ar1[0] != "-" and ar1[-1] == "-" and ar2[0] == "-" and ar2[-1] != "-":
        left = 0
        for x in ar2:
            if x == "-":
                left += 1
            else:
                break

        right = 0
        for x in ar1[::-1]:
            if x == "-":
                right += 1
            else:
                break

        length = len(ar1) - left - right

        ar1 = ar1[left:left + length]
        ar2 = ar2[left:left + length]
    else:
        # print(ar1)
        # print(ar2)
        raise WeirdReads

    # Create a canonical read.
    canonical = ""
    pairMismatches = 0

    for i in range(len(ar1)):
        nt_1 = ar1[i]
        nt_2 = ar2[i]

        # If both reads have the same nucleotide - use it
        if nt_1 == nt_2:
            canonical += nt_1
        else:
            # No match means mismatch (large insertions/deletions are counted multiple times)
            pairMismatches += 1

            # Gap on nt_1 means nt_2 must have something there.
            if nt_1 == "-":
                canonical += nt_2
            # Gap on nt_2 means nt_1 must have something
            elif nt_2 == "-":
                canonical += nt_1
            # If neither is a gap it means a mismatch between both. For now, we use "N" in this case.
            else:
                if prefer_mate:
                    canonical += nt_2
                else:
                    canonical += nt_1

    if pairMismatches > max_pair_mismatches:
        raise TooManyMismatchesBetweenPair

    # Prepare the result object
    result = AlignedReadResult(read_length)
    result.pairMismatches = pairMismatches

    # Align canonical read to reference
    pair_b = pairwise2.align.localcs(reference, canonical, _match_score, -10, -1, one_alignment_only=True)[0]

    a_ref = pair_b[0]
    a_can = pair_b[1]

    # Count mutations!
    r_pos = 0
    insertion_open = False
    deletion_open = False
    refMismatches = 0
    readMismatches = np.zeros(read_length)
    readMutationTarget = np.zeros((read_length, 5))
    indelHappened = 0
    insLength = 0
    maxInsLength = 0
    delLength = 0
    maxDelLength = 0

    for i in range(len(a_ref)):
        nt_ref = a_ref[i]
        nt_can = a_can[i]

        if insLength > maxInsLength:
            maxInsLength = insLength

        if delLength > maxDelLength:
            maxDelLength = delLength

        # Just continue if reference contains "N"
        if nt_ref == "N":
            r_pos += 1
            insertion_open = False
            deletion_open = False
            # If both agree we can increase r_pos and continue (= no mutation)
        elif nt_ref == nt_can:
            r_pos += 1
            insertion_open = False
            deletion_open = False
        else:
            # A gap on the reference is an insertion. Do not progress r_pos in this case.
            if nt_ref == "-":
                # Only count as one insertion
                if insertion_open is False:
                    result.indels[r_pos, 0] += 1
                    result.insertionTarget[r_pos, dna.base2Index(nt_can)] += 1
                    refMismatches += 1
                    indelHappened += 1

                    # Debug
                    #if r_pos == 49:
                    #    print(a_ref, a_can)
                insertion_open = True
                deletion_open = False
                insLength += 1
                delLength = 0
            # A gap on the canonical is a deletion. Do progress r_pos in this case
            elif nt_can == "-":
                if deletion_open is False:
                    result.indels[r_pos, 1] += 1
                    refMismatches += 1
                    indelHappened += 1
                r_pos += 1
                insertion_open = False
                deletion_open = True
                delLength += 1
                insLength = 0
            # Neither is a point mutation. progress r_pos in this case.
            else:
                result.mismatches[r_pos] += 1
                result.mutationTarget[r_pos, dna.base2Index(nt_can)] += 1

                r_pos += 1
                insertion_open = False
                deletion_open = False
                insLength = 0
                delLength = 0

        # If we reached the end of the read we can stop.
        # Any additional bases on the canonical strand are not counted.
        if r_pos == read_length:
            break

    if result.mismatchCount > max_point_mutations:
        raise TooManyMismatchesToReference

    if maxDelLength > 2:
        #print("Too long deletions (>2).")
        raise TooManyMismatchesToReference
    if maxInsLength > 2:
        #print("Too long insertions (>2).")
        raise TooManyMismatchesToReference

    # Split on codon if given
    codon = None
    if split_on_codon > 0:
        if split_on_codon < len(coordinates):
            try:
                x, y = coordinates[split_on_codon]
                codonRead = a_can[x:y]

                if codonRead in codons:
                    codon = codonRead
            except TypeError:
                print(type(x), x, type(y), y)
        else:
            print("Codon number is beyond limits.")

    result.codon = codon

    # Save other codons, too
    codons = []
    for x, y in coordinates:
        codons.append(a_can[x:y])
    result.codons = codons

    return result
