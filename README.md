# Sequence simulator for homolog_distance

This repository contains Python code for a program that generates synthetic data
for testing the [homolog_distance](https://github.com/actapia/homolog_distance)
program, which computes a genetic distance matrix for a set of individual
organisms using RNA-seq data.

With its most basic configuration, this software randomly generates a
phylogenetic tree and simulates the evolution of a specified number of randomly
generated transcripts on the tree. The output of the program is the generated
phylogenetic tree (represented in Newick format) and the generated transcripts
for each of the extant taxa in the tree.

This program is designed to be configurable via a combination of command-line
arguments and a YAML configuration file&mdash;basic usage of the program should
not require familiarity with Python or the source code of this program. However,
this sequence simulator is also designed to be extensible through custom
functions and simulation classes. For now, please refer to the
comments/docstrings in the source code to understand how to extend the
capabilities of this program to suit your needs.

## Requirements

This program has been tested with the following configuration:

* CPython 3.11.3
* Biopython 1.80
* Dendropy 4.6.1
* Joblib 1.3.2
* tqdm 4.64.1
* more_itertools 9.1.0
* PyYAML 6.0.1
* SciPy 1.10.1
* NumPy 1.25.1

Other similar configurations may work but are untested. Much older versions of
Joblib and Dendropy likely will *not* work&mdash;this software uses new features
of both as of this writing.

## Basic usage

The sequence simulator can be used without writing any Python code. A basic
configuration file is shown below. It may also be found at
[sim_config.yaml](./sim_config.yaml).

```yaml
# Extant taxa in the final tree.
taxa: 16
# Number of transcripts per taxon.
count: 50000
0# Random seed (for reproducibility)
seed: 487
# The function to use for generating the tree
tree_generator:
  name: birth_death_tree
  params:
    birth_rate: 1.0
    death_rate: 0.5
# The function to use for determining the lengths of transcripts.
seq_len_distribution:
  name: binomial
  params:
    p: 0.1
    n: 1000
    loc: 1950
# The function to use for simulating sequence evolution on the tree.
char_generator:
  name: hky85
  params:
    mutation_rate: 0.01
# The lowest log level to display.
log_level: DEBUG
# When specified, the simulation saves the generated tree to a file with this
# name.
save_tree: "phylogeny.tree"
# Template for naming transcripts.
seqid_template: "NODE_cov_{cov}_g{gene}_i{iso}"
# Whether to put each taxon's output in a separate directory.
separate_dirs: True
# The function to use for assigning coverage values to transcripts.
coverage_distribution:
  name: uniform
  params:
    min: 0
    max: 10000
```

This configuration file tells the simulator to generate a tree with 16 extant
taxa and generate and evolve 50000 transcripts. The random seed ensures
reproducibility. The simulation is configured to use the `birth_death_tree`
model for generating the tree, a `binomial` distribution for determining
transcript lengths, and a `uniform` distribution for assigning coverage values
to transcripts. 

The simulator is configured to use the `hky85` model for evolving sequences on
the tree. Messages at or above the log level `DEBUG` will be shown. The
generated tree will be saved to a file named `phylogeny.tree`. The template used
for naming the transcripts is `"NODE_cov_{cov}_g{gene}_i{iso}"`. Each taxon's
transcripts will be output to a FASTA file in a separate directory.

If the simulation configuration file is located at `sim_config.yaml`, then the
simulator can be run with the specified configuration with a command like the
one below.

```bash
python generate_simulated_data.py -c sim_config.yaml -j 0 -O sim_example
```

The command-line options `-j 0` and `-O sim_example` specify that the simulator
should be run with the maximum number of parallel jobs and should write the
output to a directory named `sim_example`, respectively.

## Simulation configuration files

The configuration for a simulation may be specified in a YAML file. The YAML
format is designed to be easily readable by both humans and software like this
simulator.

The subsections below describe the options understood by the default simulator.

### `taxa`

An non-negative integer specifying the number of extant taxa ("tips") present in
the final tree.

### `count`

A non-negative integer specifying the number of transcripts to generate per
taxon.

### `seed`

A non-negative integer used to initialize the pseudo-random number generators
for the simulation. Specifying the `seed` is useful for ensuring
reproduciblity&mdash;for a fixed seed, system, and version of the program, the
simulator should produce the same output each time the program is run.

If the `seed` is not specifying, an "unpredictable" seed will be used
automatically.

### `tree_generator`

A name and parameters for the model to use for generating the phylogenetic
tree. The table below summarizes the valid tree models understood by the basic
simulator.

| Model name                  | Description                                                                     | Parameters                                                                                                       |
|-----------------------------|---------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------|
| `birth_death_tree`          | A simple model in which births and deaths can occur randomly at each time step. | See [Dendropy documentation](https://dendropy.org/library/birthdeath#dendropy.model.birthdeath.birth_death_tree) |
| `extended_birth_death_tree` | Like the birth-death model, but the tree is extended after initial generation.  | See [`tree_gen.py`](tree_gen.py)                                                                                 |

### `seq_len_distribution`

A name and parameters for the function to use for randomly determining
transcript sequence lengths.

| Function name | Description                                                  | Parameters                                                                                             |
|---------------|--------------------------------------------------------------|--------------------------------------------------------------------------------------------------------|
| `binomial`    | A binomial distribution.                                     | See [SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom.html) |
| `match`       | Match the distribution of sequence lengths for the fit data. (See [`fit_data`](###fit_data).) | None.                                                                                                  |

### `char_generator`

A name and parameters for the model to use for simulating evolution of
transcripts along a tree.

| Model name | Description                                          | Parameters                                                                                              |
|------------|------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
| `hky85`    | Basic model of sequence evolution by Hasegawa et al. | See [Dendropy documentation](https://dendropy.org/library/discrete#dendropy.model.discrete.hky85_chars) |

### `log_level`

The minimum log level for which messages will be shown.

| Log level  | Description                                                             |
|------------|-------------------------------------------------------------------------|
| `DEBUG`    | Messages mainly useful for debugging the program.                       |
| `INFO`     | Useful information that does not indicate a problem.                    |
| `WARNING`  | Information related to possible errors/problems.                        |
| `ERROR`    | Information related to definite errors/problems.                        |
| `CRITICAL` | Information related to errors that prevent the program from proceeding. |

### `save_tree`

A string specifying the name of the file to which the generated phylogenetic
tree will be saved.

If no value for this option is specified, the program will *not* save the
phylogenetic tree to a file.

### `seqid_template`

A Python template string that will be used to name the transcripts generated.

The basic simulation program recognizes three fields that will be replaced with
appropriate values when their names are present in the template string
surrounded by curly braces.

| Field name | Description                       |
|------------|-----------------------------------|
| `cov`      | Coverage value of the transcript. |
| `gene`     | Gene ID.                          |
| `iso`      | Isoform ID.                       |

If no value for this option is specified, the default is `"transcript_{gene}"`.

### `separate_dirs`

A boolean indicating whether the transcripts for each taxon should be in
separate directories.

If no value for this option is specified, the default is `False`.

### `coverage_distribution`

A name and parameters for the function to use for randomly determining
transcript coverage values.

| Function name | Description             | Parameters                                                   |
|---------------|-------------------------|--------------------------------------------------------------|
| `uniform`     | A uniform distribution. | See [`coverage_distributions.py`](coverage_distributions.py) |

### `fit_data`

A string or list of strings specifying paths to data for which the `taxa` and
`seq_len_distribution` settings will be emulated.

For example, if there are ten taxa in the fit data, then the simulator will
generate ten taxa (unless `taxa` is also specified manually in the
configuration).

The `fit_data` only affects the sequence length distribution when
`seq_len_distribution` is set to `match`. In that case, the lengths of the
transcripts in the input data are counted, and lengths are selected at random
from the seen transcript lengths. The probability that any length is selected is
proportional to the number of times it appears as the length of a transcript in
the input data.

### `jobs`

A non-negative integer indicating the number of parallel jobs to use for the
simulation.

Some parts of the simulation can be parallelized to improve speed. This
parameter lets you select how many tasks will be run in parallel in such parts
of the simulation.

If the special value `0` is given for this setting, the number of jobs used will
be the number of processor cores detected on the computer running the
simulation.

### `out_dir`

A string specifying a path to the directory in which the data should be saved.


## Command-line arguments

Some options can additionally be specified via command-line arguments to the
simulation script. The table below describes the accepted command-line arguments
and the configuration file options to which they correspond. 

| Short name | Long name     | Corresponding option           |
|------------|---------------|--------------------------------|
| `-O`       | `--out-dir`   | [`out_dir`](###`out_dir`)     |
| `-f`       | `--fit-data`  | [`fit_data`](###`fit_data`)   |
| `-n`       | `--count`     | [`count`](###`count`)         |
| `-j`       | `--jobs`      | [`jobs`](###`jobs`)           |
| `-t`       | `--taxa`      | [`taxa`](###`taxa`)           |
|            | `--log-level` | [`log_level`](###`log_level`) |
|            | `--seed`      | [`seed`](###`seed`)           |

Additionally, the following two arguments do not correspond to any configuration
file options.

| Short name | Long name  | Description                      |
|------------|------------|----------------------------------|
| `-c`       | `--config` | Path to YAML configuration file. |
| `-h`       | `--help`   | Display a help message.          |

The `--help` option may be used to print a summary of these accepted
command-line arguments, but the help does not currently explain the accepted
configuration file options.
