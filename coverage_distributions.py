import random

def uniform(rng: random.Random, **kwargs) -> float:
    """Select a floating-point number uniformly in the specified range.

    Although this function accepts a variable number of arbitrary keyword
    arguments, only two can actually be used: min and max.

    Parameters:
        rng: A Python pseudo-random number generator.
        min: The minimum value that may be generated, inclusive.
        max: The maximum value that may be generated, inclusive.

    Returns:
        A float choosen uniformly at random from the range [min, max].
    """
    return rng.uniform(kwargs["min"], kwargs["max"])
