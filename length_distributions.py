import functools
import random
from collections import Counter
from numbers import Number
from typing import Iterable, TextIO, Union, Callable, TypeVar

from Bio import SeqIO

T = TypeVar("T")

def weighted_choice(weight_dict: dict[T, Number], rng: random.Random) -> T:
    """Select at random from weight_dict's keys using the values as weights

    Parameters:
        weight_dict (dict): A dictionary mapping choices to their weights.
        rng:                A Python RNG object to use for making the choice.

    Returns:
        A key selected at random from weight_dict with its value as weight.
    """
    # Thanks to pbsds on Stack Overflow for this trick.
    # https://stackoverflow.com/a/54494915
    return rng.choices(*zip(*weight_dict.items()))[0]

def match_dist(fit_data: Iterable[Union[str, TextIO]]) -> Callable[[], int]:
    """Return the length distribution of the sequences in the FASTA files.

    This function counts the distinct lengths of all sequences in all provided
    input FASTA files and returns a function that selects from those lengths at
    random in such a way that a length is chosen with probability proportional
    to the number of sequences encountered with that length.

    Parameters:
        fit_data: A list of paths or file-like objects for FASTA files.

    Returns:
        A function that matches the length distribution of the input sequences.
    """    
    return functools.partial(
        weighted_choice,
        Counter(
            len(s) for f in fit_data for s in SeqIO.parse(f, format="fasta")
        )
    )

def constant(x: T) -> Callable[..., T]:
    """Return a variadic function that returns a specified constant value."""
    def inner(*args, **kwargs):
        return x
    return inner
