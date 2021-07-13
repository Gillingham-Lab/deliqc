import argh
from colorama import init, Fore, Style
import os
from deliqc.cli.helpers import critical, warning, expand_filepattern
from deliqc.sample.load import load_sample
import glob
from deliqc import dna
from multiprocessing import cpu_count
from timeit import default_timer as timer
import numpy as np
import pickle


@argh.arg("--threads", type=int, default=4)
@argh.arg("--max-reads", type=int)
@argh.arg("--max-pair-mismatches", type=int)
@argh.arg("--max-point-mutations", type=int)
@argh.arg("--split-on-codon", type=int)
@argh.arg("--codons", type=list, nargs="+", default=[])
@argh.arg("--title", type=str)
def run(
        sequence: "DNA template sequence. Use N to mark codons.",
        r1: "Main read file to run deliqc on. Reads must be in same order as given in sequence. Use * for wildcards to select multiple files.",
        r2: "Mate file for paired end sequencing. Use * for wildcards to select multiple files.",
        save_as: "Target filename to save the pickle with extracted data for post-analysis",
        threads: "Number of threads" = None,
        max_reads: "Maximum amount of reads" = 10_000,
        max_pair_mismatches: "Maximum amount of mismatches between pairs" = 0,
        max_point_mutations: "Maximum amount of point mutations before a read is skipped" = 5,
        split_on_codon: "Codon to split the results on" = 0,
        codons: "Codons to split the results on" = [],
        title: "A title to use within the file" = None,
):
    """ Main runner for extraction"""
    init()

    print("sequence", sequence)
    print("r1", r1)
    print("r2", r2)

    # Sanitise DNA sequence and check if its valid
    sanitised_sequence = dna.sanitise(sequence)

    if sequence != sanitised_sequence:
        warning(f"Sequence was sanitised: {sanitised_sequence}")

    if not dna.is_valid(sanitised_sequence):
        critical("DNA must only contain A, T, G, C or N.")

    # Check if r1 and r2 make sense
    try:
        r1 = expand_filepattern(r1)
        r2 = expand_filepattern(r2)
    except FileNotFoundError as e:
        critical(e)
        return

    if len(r1) != len(r2):
        critical(f"r1 and r2 must have an equal amount of files.")

    # Check other parameters
    if max_reads <= 0:
        critical(f"Maximum amount of reads must be at least 1.")

    if threads <= 0:
        warning("Number of threads was less than 1 and therefore forced to 1.")
        threads = 1
    elif threads > cpu_count():
        warning("Number of threads exceeds cpu count.")

    start = timer()

    # Iterate over the samples
    print(f"Found {len(r1)} samples.")
    results = []
    for f1, f2 in zip(r1, r2):
        print(f" - Starting to extract {max_reads} from {f1} and {f2}")
        result = load_sample(
            sanitised_sequence,
            f1, f2,
            max_reads=max_reads,
            max_pair_mismatches=max_pair_mismatches,
            max_point_mutations=max_point_mutations,
            split_on_codon=split_on_codon,
            codons=codons
        )

        print(f"    - Reads found: {result['reads']}")
        print(f"    - Reads aligned: {result['alignedReads']}")
        print(f"    - Reads skipped (mate mismatches): {result['mismatchedPairs']}")
        print(f"    - Reads skipped (too many point mutations): {result['tooManyPointMutations']}")
        print()

        results.append(result)

    sample = {
        "title": title if title is not None else r1[0].split(".")[0],
        "reference": sanitised_sequence,
        "splitOnCodon": split_on_codon,
        "splitOnCodonSequences": codons,
        "replicates": results,
        "mismatches": (
            np.mean([r["mismatches"] for r in results], axis=0),
            np.std([r["mismatches"] for r in results], axis=0),
        ),
        "indels": (
            np.mean([r["indels"] for r in results], axis=0),
            np.std([r["indels"] for r in results], axis=0),
        ),
        "mutationTarget": (
            np.mean([r["mutationTarget"] for r in results], axis=0),
            np.std([r["mutationTarget"] for r in results], axis=0),
        ),
        "insertionTarget": (
            np.mean([r["insertionTarget"] for r in results], axis=0),
            np.std([r["insertionTarget"] for r in results], axis=0),
        ),

        "reads": np.mean([r["reads"] for r in results]),
        "mismatchedPairs": np.mean([r["mismatchedPairs"] for r in results]),
        "tooManyPointMutations": np.mean([r["tooManyPointMutations"] for r in results]),
        "alignedReads": np.mean([r["alignedReads"] for r in results]),
    }

    aligned_reads_per_codon = {}
    for k in results[0]["alignedReadsPerCodon"].keys():
        aligned_reads_per_codon[k] = np.mean([r["alignedReadsPerCodon"][k] for r in results])
    sample["alignedReadsPerCodon"] = aligned_reads_per_codon

    print(f" - Average reads found: {sample['reads']:.0f}")
    print(f" - Average reads aligned: {sample['alignedReads']:.0f}")
    print(f" - Average reads skipped (mate mismatches): {sample['mismatchedPairs']:.0f}")
    print(f" - Average reads skipped (too many point mutations): {sample['tooManyPointMutations']:.0f}")
    print(f" - Average point (potential priming site): {sample['mismatches'][0][:8, 0].mean()*100:.1f}%")
    print(f" - Average point mutation rate: {sample['mismatches'][0][:, 0].mean()*100:.1f}%")
    print(f" - Average point insertion rate: {sample['indels'][0][:, 0, 0].mean()*100:.1f}%")
    print(f" - Average point deletion rate: {sample['indels'][0][:, 1, 0].mean()*100:.1f}%")

    stop = timer()

    print(f"\n... Done [{stop-start:.1f} seconds]")

    with open(save_as, "wb") as fh:
        pickle.dump(sample, fh)

    print(f"File saved ({save_as})")
