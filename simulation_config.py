import argparse
import functools
import glob
import logging
import multiprocessing
import operator
from pathlib import Path
from typing import Callable, Union, Any, TextIO

import yaml
from Bio import SeqIO
from dendropy.simulate import treesim
from dendropy.model import discrete
from scipy import stats

from tree_gen import extended_birth_death_tree
from length_distributions import match_dist
from coverage_distributions import uniform
from simulation import Simulation

def int_or_str(valid: set[str]) -> Callable[[Any], Union[int, str]]:
    """Returns a function that parses its argument as an int or one valid str.

    This function is useful for parsing arguments that can either take an
    integer value or one of a set number of valid strings.

    Parameters:
        valid (set): A set of valid strings to allow.

    Returns:
        A function that parses its parameter as an integer or a string in valid.
    """    
    def inner(v):
        try:
            return int(v)
        except ValueError as e:
            if v in valid:
                return v
            else:
                raise e
    return inner

def count_sequences(fit_data: Union[str, TextIO]) -> int:
    """Returns the number of sequences in the given FASTA file."""
    return max(sum(1 for _ in SeqIO.parse(f, format="fasta")) for f in fit_data)

log_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]

def use_np(f: Callable, argname: str = "random_state") -> Callable:
    """Returns a function that passes the NumPy RNG to the given function.

    This function and use_py are useful to allow functions requiring a
    pseudo-random number generator be called using the same parameters
    regardless of whether they require a NumPY RNG or a Python RNG.

    This function accepts a function f as its first parameter. This function
    returns a function inner that calls f with a keyword argument specified by
    this function's argname parameter set to the NumPy RNG passed to inner.

    Parameters:
        f:             The function to pass a NumPy RNG.
        argname (str): The formal parameter for the NumPy RNG accepted by f.

    Returns:
        A function that takes Python and NumPy RNGs and passes the second to f.
    """
    def inner(*args, py_rand, np_rand, **kwargs):
        kwargs[argname] = np_rand
        return f(*args, **kwargs)
    return inner

def use_py(f: Callable, argname: str = "rng") -> Callable:
    """Returns a function that passes the NumPy RNG to the given function.

    This function and use_np are useful to allow functions requiring a
    pseudo-random number generator be called using the same parameters
    regardless of whether they require a NumPy RNG or a Python RNG.

    This function accepts a function f as its first parameter. This function
    returns a function inner that calls f with a keyword argument specified by
    this function's argname parameter set to the Python RNG passed to inner.

    Parameters:
        f:             The function to pass a Python RNG.
        argname (str): The formal parameter for the Python RNG accepted by f.

    Returns:
        A function that takes Python and NumPy RNGs and passes the first to f.
    """    
    def inner(*args, py_rand, np_rand, **kwargs):
        kwargs[argname] = py_rand
        return f(*args, **kwargs)
    return inner

def preset(f: Callable) -> Callable[[dict], Callable]:
    """Returns a function that creates a partially-applied version of f.

    This function accepts a function f as its only parameter and returns a
    function inner that accepts a dictionary of keyword arguments used to
    create and return a partially applied version of f.

    This function is similar to joblib's delayed function, but the function it
    returns accepts a dictionary of keyword arguments rather than a variable
    number of arguments and keyword arguments.

    Parameters:
        f: A function.

    Returns:
        A function that returns f with the provided arguments partially applied.
    """
    def inner(kwargs):
        return functools.partial(f, **kwargs)
    return inner

def load_from_specs(d: dict[str, Callable[[dict], Any]], c: dict[str, Any]):
    """Returns the result of calling a function in d with specified parameter.

    The first parameter d must map names to functions. The second parameter c
    must specify the name of the function and the dict parameter to give to the
    function.

    Parameters:
        d: A dictionary mapping function names to functions.
        c: A dictionary containing the name and dict parameter for a function.

    Returns:
        The result of calling the function named in c with the parameter in c.
    """
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
    programs. If you simply need to create a Simulation in your own code, you
    may find it easier to simply construct a Simulation.

    Attributes:
        default_config (dict):   The default base configuration.
        config (dict):           The current configuration.
        simulation_class (type): The simulation class to construct.
    """

    # This dict specifies config options to be treated as names of functions.
    # The values of the dict are dicts mapping function names to "presets" of
    # functions that return partially applied functions.
    load_functions = {
        # Functions for generating phylogenetic trees.
        "tree_generator": {
            "birth_death_tree": preset(use_py(treesim.birth_death_tree)),
            "extended_birth_death_tree": preset(
                use_py(extended_birth_death_tree)
            ),
        },
        # Functions for simulating sequence evolution on trees.
        "char_generator": {
            "hky85": preset(use_py(discrete.hky85_chars)),
        },
        # Functions for selecting sequence lengths for transcripts.
        "seq_len_distribution": {
            "binomial": lambda kwargs: use_np(stats.binom(**kwargs).rvs),
            "match": lambda kwargs: use_py(match_dist(**kwargs))
        },
        # Functions for selecting coverage values for transcripts.
        "coverage_distribution": {
            "uniform": preset(use_py(uniform))
        }
    }
    
    def __init__(
            self,
            default_config: dict[str, Any] = None,
            simulation_class: type = Simulation
    ):
        """Construct a config with specified defaults and simulation type.

        Parameters:
            default_config (dict):   The base configuration to start with.
            simulation_class (type): The simulation class to construct.
        """
        if default_config is None:
            default_config = {}
        self.config = {}
        self.default_config = default_config
        self.simulation_class = simulation_class
        
    def _make_parser(self) -> argparse.ArgumentParser:
        """Create a command-line argument parser for this configuration."""
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

    def _handle_arguments(self) -> argparse.Namespace:
        """Parse arguments from the command line."""
        return self._make_parser().parse_args()

    def _set_logging(self):
        """Set the logging level to the configured value."""
        logging.basicConfig(level=getattr(logging, self.config["log_level"]))

    def _pre_cli(self):
        """Handle configuration settings before parsing command-line args."""
        if "fit_data" in self.config:
            try:
                self.config["fit_data"] = glob.glob(self.config["fit_data"])
            except TypeError:
                self.config["fit_data"] = functools.reduce(
                    operator.add,
                    (glob.glob(f) for f in self.config["fit_data"])
                )

    def _cli(self, args: argparse.Namespace):
        """Add command-line arguments to the configuration."""
        args_dict = {
            k: v for (k, v) in vars(args).items() if v is not None
        }
        cli_options = set(args_dict.keys()) - {"config"}
        for option in cli_options:
            self.config[option] = args_dict[option]

    def _handle_load_functions(self):
        """Load functions named in configuration as Python functions.

        Since YAML does not provide an easy way to set a configuration key to
        a Python function, this class interprets certain pre-defined keys as
        referring to functions.

        A value for one of these keys is itself a dictionary with the key
        "name", used to identify which function is associated with the first
        key, and "params", which is a dictionary of keyword arguments to be
        passed to the function identified by the "name" key.

        See the load_from_specs function and the load_functions attribute of
        this class to understand in more detail how these functions are loaded
        according to the options specified in a config file.
        """
        for f, d in self.load_functions.items():
            try:
                self.config[f] = load_from_specs(d, self.config[f])
            except KeyError:
                pass

    def _post_cli(self):
        """Handle configuration settings after parsing command-line args."""
        if self.config["count"] == "match":
            logging.debug("Counting sequences.")
            self.config["count"] = count_sequences(self.config["fit_data"])
        if self.config["taxa"] == "match":
            self.config["taxa"] = len(self.config["fit_data"])
        # if "taxa" in self.config:
        #     if self.config["tree_generator"]["name"] == "birth_death_tree":
        #         self.config["tree_generator"]["params"].setdefault(
        #             "num_extant_tips",
        #             self.config["taxa"]
        #         )
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

    def load_config(
            self,
            args: argparse.Namespace = None,
            config: Union[str, Path] = None,
            set_logging: bool = True
    ):
        """Load the configuration from the given sources.

        This function loads a configuration from up to two
        sources---command-line arguments and a configuration file.

        Parameters:
            args:               An Namespace containing parsed CLI arguments.
            config:             A path to a YAML configuration file to load.
            set_logging (bool): Whether to set the logging level.
        """
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
        """Construct a simulation object using this configuration."""
        sim = self.simulation_class()
        for k, v in self.config.items():
            setattr(sim, k, v)
        sim.config = self
        return sim
