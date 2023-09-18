# Hybridization and divergence simulation

This repository also includes code for a "hybridization and divergence"
simulation; this code also serves as an example of extending the base simulation
code for a custom simulation.

The hybridization and divergence simulation involves the following steps:

1. Generate a phylogenetic tree and simulate evolution of randomly generated
sequences on the tree. (This is the same as the base simulation.)
	   
2. Selection of random disjoint taxon sets of specified size to "hybridize."
Transcripts of the taxa in the set are combined into a single file,
renaming the transcripts to identify their origin. This produces a new
"hybrid" taxon.
	   
3. Further evolution of the hybrid taxa by generating a new tree for each
hybrid and simulating evolution of the transcripts along that new tree.
	   
The intermediate results of the first two steps are saved to their own
directories. New configuration options may be used to specify the number of taxa
to include in a hybrid, the number of descendant taxa in the new trees, and the
output directories.

# Basic usage

The example configuration file below can be used with the hybridization and
divergence simulation. You can also find this file at
[divergence_example.yaml](../divergence_example.yaml).

```yaml
# Basic simulation params.
# See the main README.md for explanations of these parameters.

taxa: 16
count: 50000
seed: 487
tree_generator:
  name: extended_birth_death_tree
  params:
    birth_rate: 1.0
    death_rate: 0.5
    extend_relative: 0.1
seq_len_distribution:
  name: binomial
  params:
    p: 0.1
    n: 1000
    loc: 1950
char_generator:
  name: hky85
  params:
    mutation_rate: 0.01
log_level: DEBUG
save_tree: "phylogeny.tree"
seqid_template: "NODE_cov_{cov}_g{gene}_i{iso}"
separate_dirs: True
coverage_distribution:
  name: uniform
  params:
    min: 0
    max: 10000


# Hybrid params

# The number of taxa to combine in each hybrid.
select: 2


# Divergence params

# The number of extant taxa in the new trees to generate after hybridization.
diverged_taxa: 4
```

The configuration above includes many familiar options that are described in the
[README.md](../REAMDE.md) at the root of this repository.

Assuming the above configuration is in a file named `divergence_example.yaml`,
the hybridization and divergence simulation can be run with the following
command.

```bash
python do_hybrid_divergence_sim.py \
    -c divergence_example.yaml \
	-j 0 \
	--divergence-out-root hybrid5
```

The `-c` and `-j` command-line arguments have the same meaning as in the basic
simulation. The `--divergence-out-root` option specifies that the output
directories for all three phases of the simulation should be created in the
`hybrid5` directory.

After running the command, the `hybrid5` directory will contain three
subdirectories&mdash;`p1`, `p2`, and `p3`, for the first, second, and third
phases of the simulation, respectively.

### p1 directory

The `p1` directory contains the outputs of the first part of the simulation. Its
structure is identical to that of the output from the base simulation. If the
simulation is configured to save the phylogenetic trees generated, a tree with
the configured filename can be found in the root of the `p1` directory.

The `p1` directory also contains the transcripts for the extant taxa of the
generated tree. If the `separate_dirs` option was set to `True` in the
configuration, then each transcripts FASTA file will reside in a directory
specific to the taxon to which the transcripts belong.

### p2 directory

The `p2` directory contains the output of the second part of the
simulation&mdash;i.e., the transcripts for the hybrid taxa. The hybrid taxa and
their transcripts are identified by their parent taxa. If `separate_dirs` is on,
then each directory/taxon has a name consisting of one or more taxon names from
the first simulation phase separated by a plus sign. 

### p3 directory

The `p3` directory contains subdirectories for each of the hybrid taxa generated
in phase 2. These subdirectories contain the transcripts of the descendants of
the hybrids. Again, if `save_tree` is specified, there will be a phylogenetic
tree with the specified name in each of these subdirectories.

The names of the descendant taxa consist of a prefix denoting the parent
ancestral hybrid followed by a descendant ID beginning with "T". Again, if
`separate_dirs` is on, then each descendant taxon will get its own subdirectory.

## Configuration options

In addition to the configuration options supported by the base simulation, the
hybridization and divergence simulation supports the options described below.

(For information on the base simulation options, refer to the
[README.md](../README.md) at the root of this repository.)

### `select`

A positive integer representing the number of taxa to select for creating each
hybrid.

Since taxa are selected for hybridization without replacement, `taxa` must be
divisible by `select`.

### `hybrid_out_dir`

A string specifying where to put the output files from the second phase of the
simulation.

If the [`divergence_out_root`](#divergence_out_root) option is specified, the
`hybrid_out_dir` will automatically be set to a subdirectory `p2` of the
specified `divergence_out_root` unless `hybrid_out_dir` is manually specified.

### `diverged_taxa`

A non-negative integer representing the number of new extant taxa in the
generated trees in the last phase.

Currently, there is no way to have a different number of extant taxa for each
new tree.

### `divergence_out_dir`

A string specifying where to put the output files from the third phase of the
simulation.

If the [`divergence_out_root`](#divergence_out_root) option is specified, the
`divergence_out_dir` will automatically be set to a subdirectory `p3` of the
specified `divergence_out_root` unless `divergence_out_dir` is manually
specified.

### `divergence_out_root`

A string specifying the directory in which to put the output directories for the
three simulation phases.

If this option is specified, then the options
[`out_dir`](../README.md#out_dir), [`hybrid_out_dir`](#hybrid_out_dir), and
[`divergence_out_dir`](#divergence_out_dir) will be set automatically to
subdirectories `p1`, `p2`, and `p3` of the `divergence_out_root`, respectively.

## Command-line arguments

The new command-line arguments specified in the table below are recognized. All
correspond to options that may be specified in a configuration file.

| Short name | Long name               | Corresponding option                          |
|------------|-------------------------|-----------------------------------------------|
|            | `--hybrid-out-dir`      | [`hybrid_out_dir`](#hybrid_out_dir)           |
|            | `--divergence-out-dir`  | [`divergence_out_dir`](#divergence_out_dir)   |
|            | `--divergence-out-root` | [`divergence_out_root`](#divergence_out_root) |

As in the default simulation, you can specify the `--help` option to see an
overview of the command-line arguments accepted by the
`do_hybrid_divergence_sim.py` script.
