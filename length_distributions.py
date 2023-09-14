import functools
from collections import Counter

from Bio import SeqIO

def weighted_choice(weight_dict, rng):
    # Thanks to pbsds on Stack Overflow for this trick.
    # https://stackoverflow.com/a/54494915
    return rng.choices(*zip(*weight_dict.items()))[0]

def match_dist(fit_data):
    return functools.partial(
        weighted_choice,
        Counter(
            len(s) for f in fit_data for s in SeqIO.parse(f, format="fasta")
        )
    )
