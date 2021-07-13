import gzip


def read_lines(source: str, max_reads: int =20):
    """
    Reads a given amount of lines from a source.

    Source can be either a fastq file (.fq) or a gzip compressed fastq file (.fq.gz or .gz).
    :param source:
    :param max_reads:
    :return:
    """
    if source.endswith("gz"):
        h = gzip.open
        decode = True
    else:
        h = open
        decode = False

    with h(source, "r") as fh:
        worked_lines = 0
        for header in fh:
            read = next(fh)
            qs = next(fh)
            qual = next(fh)

            if decode:
                yield read.strip().decode("ascii")
            else:
                yield read.strip()

            worked_lines += 1

            if worked_lines >= max_reads:
                break
