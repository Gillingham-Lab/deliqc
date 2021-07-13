from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, PercentFormatter, FuncFormatter)
import numpy as np

from deliqc import dna
from .rcparams import get_rcparams


def plot_point_mutation_rate(sample, codon_index=0, figsize=(6, 1.5), dpi=300):
    ci = codon_index

    title = sample["title"]
    reference = sample["reference"]
    L = len(reference)

    x = range(1, L + 1)
    y = sample["mismatches"][0][:, ci]
    y_err = sample["mismatches"][1][:, ci]
    color_sequence = dna.get_nucleotide_colour_sequence(reference)

    with plt.style.context(get_rcparams()):
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

        ax.bar(x, y*100, color=color_sequence, yerr=y_err)

        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(1))

        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.0f} %"))

        ax.set_xlim(0, L+1)
        ax.set_ylim(0, np.ceil(y.max()*100))

        ax.set_xlabel("position")
        ax.set_ylabel("point mutation rate")

        plt.tight_layout()
        plt.show()
