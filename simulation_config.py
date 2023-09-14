import argparse
import functools
import glob
import logging
import multiprocessing
import operator
from pathlib import Path

import yaml
from Bio import SeqIO
from dendropy.simulate import treesim
from dendropy.model import discrete
from scipy import stats

from tree_gen import extended_birth_death_tree
from length_distributions import match_dist
from coverage_distributions import uniform
from simulation import Simulation

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

def count_sequences(fit_data):
    return max(sum(1 for _ in SeqIO.parse(f, format="fasta")) for f in fit_data)

log_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]

def use_np(f, argname="random_state"):
    def inner(*args, py_rand, np_rand, **kwargs):
        kwargs[argname] = np_rand
        return f(*args, **kwargs)
    return inner

def use_py(f, argname="rng"):
    def inner(*args, py_rand, np_rand, **kwargs):
        kwargs[argname] = py_rand
        return f(*args, **kwargs)
    return inner

def preset(f):
    def inner(kwargs):
        return functools.partial(f, **kwargs)
    return inner

def load_from_specs(d, c):
    return d[c["name"]](c.get("params", {}))



class SimulationConfig:
    """A configuration manager for creating Simulation objects.

    This class allows a "hybrid" approach to configuring a Simulation according
    to user specifications---configuration may be performed via a config file,
    command-line arguments, or both.

    "Loading" a configuration involves the following steps:
    
        1. Considering the options that must be dealt with BEFORE the
           command-line options.

        2. Considering the command-line options.

        3. Considering the options that must be dealt with AFTER the
           command-line options.

    Each of these steps is a "private" method of this class and may be
    overridden to handle additional options in classes inheriting from this one.

    This class is mainly designed for writing the user-interfaces of simulation
    programs. If you simply need to create a Simulation in your own code, yo
    may find it easier to simply construct a Simulation.
    """

    # This dict specifies config options to be treated as names of functions.
    # The values of the dict are dicts mapping function names to "presets" of
    # functions that return partially applied functions.
    load_functions = {
        "tree_generator": {
            "birth_death_tree": preset(use_py(treesim.birth_death_tree)),
            "extended_birth_death_tree": preset(
                use_py(extended_birth_death_tree)
            ),
        },
        "char_generator": {
            "hky85": preset(use_py(discrete.hky85_chars)),
        },
        "seq_len_distribution": {
            "binomial": lambda kwargs: use_np(stats.binom(**kwargs).rvs),
            "match": lambda kwargs: use_py(match_dist(**kwargs))
        },
        "coverage_distribution": {
            "uniform": preset(use_py(uniform))
        }
    }
    
    def __init__(self, default_config=None, simulation_class=Simulation):
        if default_config is None:
            default_config = {}
        self.config = {}
        self.default_config = default_config
        self.simulation_class = simulation_class
        
    def _make_parser(self):
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
        parser.add_argument(
            "--jobs",
            "-j",
            type=int
        )
        return parser

    def _handle_arguments(self):
        return self._make_parser().parse_args()

    def _set_logging(self):
        logging.basicConfig(level=getattr(logging, self.config["log_level"]))

    def _pre_cli(self):
        if "fit_data" in self.config:
            try:
                self.config["fit_data"] = glob.glob(self.config["fit_data"])
            except TypeError:
                self.config["fit_data"] = functools.reduce(
                    operator.add,
                    (glob.glob(f) for f in self.config["fit_data"])
                )

    def _cli(self, args):
        args_dict = {
            k: v for (k, v) in vars(args).items() if v is not None
        }
        cli_options = set(args_dict.keys()) - {"config"}
        for option in cli_options:
            self.config[option] = args_dict[option]

    def _handle_load_functions(self):
        for f, d in self.load_functions.items():
            try:
                self.config[f] = load_from_specs(d, self.config[f])
            except KeyError:
                pass

    def _post_cli(self):
        if self.config["count"] == "match":
            logging.debug("Counting sequences.")
            self.config["count"] = count_sequences(self.config["fit_data"])
        if self.config["taxa"] == "match":
            self.config["taxa"] = len(self.config["fit_data"])
        if "taxa" in self.config:
            if self.config["tree_generator"]["name"] == "birth_death_tree":
                self.config["tree_generator"]["params"].setdefault(
                    "num_extant_tips",
                    self.config["taxa"]
                )
        if self.config["seq_len_distribution"]["name"] == "match":
            self.config["seq_len_distribution"].setdefault("params", {}).setdefault(
                "fit_data",
                self.config.get("fit_data")
            )
        if self.config["jobs"] <= 0:
            self.config["jobs"] = (
                (self.config["jobs"] - 1)%multiprocessing.cpu_count()
            ) + 1
        self.config.setdefault("taxa_prefix", "T")
        # Load functions.
        self._handle_load_functions()

    def load_config(self, args=None, config=None, set_logging=True):
        if args is None:
            args = self._handle_arguments()
        if config is None:
            config = args.config
        with open(config, "r") as conf_file:
            self.config = self.default_config | yaml.safe_load(conf_file)
        if set_logging:
            # Set logging level.
            self._set_logging()
        # PRE-cli changes
        self._pre_cli()
        # Add CLI arguments
        self._cli(args)
        # POST-cli changes
        self._post_cli()

    def make_simulation(self):
        sim = self.simulation_class()
        for k, v in self.config.items():
            setattr(sim, k, v)
        sim.config = self
        return sim

    # def _load_config(self):
    #     return self._get_config()
