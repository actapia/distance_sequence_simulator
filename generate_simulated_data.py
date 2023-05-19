import argparse
import operator
import random
import glob
import functools
import contextlib
import logging
from pathlib import Path
from collections import Counter

from tqdm import tqdm

import yaml

from IPython import embed

import numpy as np
from scipy import stats

from dendropy.simulate import treesim
from dendropy.model import discrete
from dendropy.utility import GLOBAL_RNG as DENDROPY_RNG

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_sequences(
        tree,
        seq_len_dist,
        generate_chars,
        count=1
):
    # lens = seq_len_dist.rvs(size=count)
    for _ in tqdm(range(count)):
        #embed()
        dna_chars = generate_chars(seq_len_dist(), tree)
        yield dna_chars

def seeded(f):
    def inner(*args, seed, **kwargs):
        seed_random(seed)
        f(*args, **kwargs)
    return inner

# New idea: Use a configuration file to control generation.
# Some parameters could also be optionally controlled via command line args.

def int_or_str(valid):
    def inner(v):
        try:
            return int(v)
        except ValueError as e:
            if v in valid:
                return v
            else:
                raise e
    return inner

log_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n",
        "--count",
        type=int_or_str({"match"}),
        help="number of transcripts to generate per taxon"
    )
    parser.add_argument(
        "-t",
        "--taxa",
        type=int_or_str({"match"}),
        help="number of taxa to generate",
    )
    parser.add_argument("-O", "--out-dir", type=Path)
    parser.add_argument("-f", "--fit-data", nargs="+", type=Path)
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        default=Path("sim_config.yaml")
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="seed to use for pseudorandom number generation"
    )
    parser.add_argument(
        "--log-level",
        choices=log_levels
    )
    return parser.parse_args()

default_config = {
    "count": "match",
    "taxa": "match",
    "log_level": "WARNING",
    "log_transcript_every": 1000,
    "seqid_template": "transcript_{gene}",
    "separate_dirs": False
}

def get_config(args, set_logging=True):
    with open(args.config, "r") as conf_file:
        config = default_config | yaml.safe_load(conf_file)
    if set_logging:
        # Set logging level.
        logging.basicConfig(level=getattr(logging, config["log_level"]))
    # PRE-cli changes
    if "fit_data" in config:
        try:
            config["fit_data"] = glob.glob(config["fit_data"])
        except TypeError:
            config["fit_data"] = functools.reduce(
                operator.add,
                (glob.glob(f) for f in config["fit_data"])
            )
    args_dict = {k: v for (k, v) in vars(args).items() if v is not None}    
    cli_options = set(args_dict.keys()) - {"config"}
    for option in cli_options:
        config[option] = args_dict[option]
    # POST-cli changes
    if config["count"] == "match":
        logging.debug("Counting sequences.")
        config["count"] = count_sequences(config["fit_data"])
    if config["taxa"] == "match":
        config["taxa"] = len(config["fit_data"])
    if "taxa" in config:
        if config["tree_generator"]["name"] == "birth_death_tree":
            config["tree_generator"]["params"].setdefault(
                "num_extant_tips",
                config["taxa"]
            )
    if config["seq_len_distribution"]["name"] == "match":
        config["seq_len_distribution"].setdefault("params", {}).setdefault(
            "fit_data",
            config.get("fit_data")
        )
    return config

def preset(f):
    def inner(kwargs):
        return functools.partial(f, **kwargs)
    return inner

def weighted_choice(weight_dict):
    # Thanks to pbsds on Stack Overflow for this trick.
    # https://stackoverflow.com/a/54494915
    return random.choices(*zip(*weight_dict.items()))[0]

def match_dist(fit_data):
    return functools.partial(
        weighted_choice,
        Counter(
            len(s) for f in fit_data for s in SeqIO.parse(f, format="fasta")
        )
    )

tree_generators = {
    "birth_death_tree": preset(treesim.birth_death_tree)
}

char_generators = {
    "hky85": preset(discrete.hky85_chars)
}

seq_len_distributions = {
    "binomial": lambda kwargs: stats.binom(**kwargs).rvs,
    "match": lambda kwargs: match_dist(**kwargs)
}

def load_from_specs(d, c):
    return d[c["name"]](c.get("params", {}))

def count_sequences(fit_data):
    return max(sum(1 for _ in SeqIO.parse(f, format="fasta")) for f in fit_data)

def seed_random(seed):
    random.seed(seed)
    np.random.seed(seed + 1)
    DENDROPY_RNG.seed(seed + 2)

def main():
    config = get_config(handle_arguments())
    # Set random seed if provided.
    if "seed" in config:
        seed_random(config["seed"])
    config["out_dir"].mkdir(exist_ok=True)
    tree_gen = load_from_specs(tree_generators, config["tree_generator"])
    logging.debug("Generating tree.")
    tree = tree_gen()
    logging.debug(tree.as_ascii_plot())
    if "save_tree" in config:
        with open(Path(config["out_dir"])/config["save_tree"], "w") as out_file:
            out_file.write(tree.as_string(schema="newick"))
    #embed()
    char_gen = load_from_specs(char_generators, config["char_generator"])
    seq_len_dist = load_from_specs(
        seq_len_distributions,
        config["seq_len_distribution"]
    )
    with contextlib.ExitStack() as stack:
        logging.debug("Opening output files for writing.")
        files = []
        for t in tree.taxon_namespace:
            if config["separate_dirs"]:
                taxon_dir = config["out_dir"] / Path(t.label)
                taxon_dir.mkdir(exist_ok=True)
                files.append(
                    stack.enter_context(
                        open(taxon_dir / "transcripts.fasta", "w")
                    )
                )
            else:
                files.append(
                    stack.enter_context(
                        open(
                            config["out_dir"] / Path(
                                "{}.fasta".format(t.label)
                            ),
                            "w"
                        )
                    )
                )
        for i, mat in enumerate(
                generate_sequences(
                    tree,
                    seq_len_dist,
                    char_gen,
                    config["count"]
                )
        ):
            # if i%config["log_transcript_every"] == 0:
            #     logging.info(f"Writing transcript {i}")
            try:
                cov = random.uniform(
                    config["random_cov"].get("min", 0),
                    config["random_cov"].get("max", 10000)
                )
            except KeyError:
                cov = float(config["count"]-i)
            seq_id = config["seqid_template"].format(
                # At some point, we may simulate coverage.
                cov=cov,
                gene=i,
                # At some point, we may allow multiple isotigs.
                iso=1
            )
            #embed()
            for f, seq in zip(files, mat.values()):
                SeqIO.write(
                    SeqRecord(Seq(str(seq)), id=seq_id, description=""),
                    f,
                    format="fasta"
                )

if __name__ == "__main__":
    main()
