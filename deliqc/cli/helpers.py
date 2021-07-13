from colorama import init, Fore, Style
from glob import glob
import os
from typing import List


class Error(Exception):
    """
    Generic error exception that can get raised to let deliqc print an error.
    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        init()

        return self.msg


def critical(msg: str) -> None:
    """
    Issues an error and raises the Error exception.
    :param msg:
    :return:
    """
    raise Error(f"{Fore.RED}! {msg}{Fore.RESET}")


def warning(msg: str) -> None:
    """
    Issues a warning
    :param msg:
    :return:
    """
    print(f"{Fore.YELLOW}- {msg}{Fore.RESET}")


def expand_filepattern(pattern: str) -> List[str]:
    """
    Checks if a given file exists, and if not, tries to expand the filename using glob.
    Returns in any case a list of at least 1 file or raises a FileNotFoundError.

    :param pattern:
    :return: List of filenames
    """
    if os.path.exists(pattern):
        return [pattern]

    files = glob(pattern)

    if len(files) == 0:
        raise FileNotFoundError(f"File {pattern} was not found.")

    return files


