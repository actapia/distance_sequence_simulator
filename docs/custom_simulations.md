# Creating custom simulations

This document aims to provide information useful for creating custom simulations
using the code in this repository.

## Design overview

Most of the code in this repository exists within descendants of two
classes&mdash;`Simulation` and `SimulationConfig`. Their purposes are described
in the subsections below.

### Simulation

As the name suggests, the `Simulation` class and its descendants represent
simulations. These classes contain the logic of a simulation but do not offer
any kind of user interface.

In the base `Simulation` class, the major steps of the simulation are separated
into "private" methods that are called from the "public" `simulate` method.
This design enables descendants of `Simulation` to override the behavior for a
single part of the simulation without needing to touch any of the other code.

In "client-code," the `Simulation` object can be configured both at and after
construction&mdash;as many parameters can be specified at the constructor as one
likes. If needed, other attributes can be set after construction to alter the
settings for the `Simulation`. Note, however, that the `simulate` method does
not check that the required attributes have all been set; it is up to the client
code to ensure that no errors occur due to missing attributes.

### SimulationConfig

Again, the name of the class explains its purpose. The `SimulationConfig`
represents a configuration for a `Simulation` and can be used to create a
`Simulation` with the specified configuration by calling the `make_simulation`
method.

The main purpose of the `SimulationConfig` class is to provide a way to load
user-specified configurations for a `Simulation`. The `SimulationConfig` works
as a kind of user interface to a `Simulation`.

The type of `Simulation` that the `SimulationConfig` creates may be specified by
providing the `simulation_class` argument to the constructor. Alternatively, one
may set the `simulation_class` attribute after construction.

Loading a `SimulationConfig` from a user-specified configuration file and
command-line arguments is accomplished in three steps.

1. Handle the options in the configuration file that must be considered *before*
   the command-line arguments are considered.
   
2. Handle the command-line arguments.

3. Handle thee options in the configuration file that must be considered *after*
   the command-line arguments.
   
These three steps correspond to the `_pre_cli`, `_cli`, and `_post_cli` methods
of the `SimulationConfig` class, respectively.

### Handling functions specified in a configuration file

The base `SimulationConfig` additionally has a function `_handle_load_functions`
that transforms functions and parameters specified by name in the user's
configuration file into partially-applied Python functions that may be used in
`Simulation` objects created by the `SimulationConfig`.

The `SimulationConfig.load_functions` class attribute is a dictionary containing
keys that specify the names of the configuration keys that are interpreted as
functions. 

The values of the keys in this dictionary are themselves dictionaries containing
keys that specify the valid values for the `name` sub-attribute in the
configuration file&mdash;these are the names of the possible functions that the
user may select.

The values in these inner dictionaries are functions that return
partially-applied versions of other functions. For example,
`preset(use_py(treesim.birth_death_tree))({"birth_rate": 0.5})`, returns a
partially applied version of `use_py(treesim.birth_death_tree)` with the
`birth_rate` keyword argument fixed as `0.5`.

#### `use_py` and `use_np`

The `use_py` and `use_np` functions provide a unified interface for calling
functions that require specifying a pseudo-random number generator, regardless
of whether they use a NumPy PRNG or a Python PRNG.

The `use_py` function returns a new function that accepts the `np_rand` and
`py_rand` keyword arguments and calls the original function with the `py_rand`
argument passed as the appropriate keyword argument of the original function
(usually this is `rng`).

Likewise, the `use_np` function returns a new function that accepts the
`np_rand` and `py_rand` keyword arguments and calls the original function with
the `np_rand` argument passed as the appropriate keyword argument of the
original function (usually `random_state`).

### Specifying command-line arguments

The `_make_parser` "private" method of the `SimulationConfig` class creates an
`argparse.ArgumentParser` object that is used by the `_handle_arguments` method
to parse command-line arguments.

Subclasses of `SimulationConfig` can add new command-line arguments by calling
the `super()._make_parser()` method and modifying the returned `ArgumentParser`
object.

## Approaches to making custom simulations

There are at least three ways to extend the code in this repository to create
custom simulations. Each has certain benefits and drawbacks.

### Modifying the original code

One straightforward (perhaps obvious) approach to making a custom `Simulation`
is to simply modify the code in this repository to make the `Simulation` behave
as you like.

#### Advantages

* Straightforward
* Very flexible

#### Disadvantages

* Requires finding the part of the code to be changed.
* Makes incorporating future updates to this software more difficult.
* Makes it difficult for readers of the code to see what was changed and how.
* Changes to `Simulation` might require changes to `SimulationConfig`.

### Providing the Simulation different functions to customize behavior

Some of the behavior of a `Simulation` instance is determined by the functions
the `Simulation` object has as attributes. Specifically, the base `Simulation`
has four function attributes designed to be configured.

| Attribute               | Description                                                               |
|-------------------------|---------------------------------------------------------------------------|
| `tree_generator`        | A function that generates a random phylogenetic tree.                     |
| `char_generator`        | A function that simulates sequence evolution on a tree.                   |
| `seq_len_distribution`  | A function that assigns a (possibly random) length for each transcript.   |
| `coverage_distribution` | A function that assigns a (possibly random) coverage for each transcript. |

By providing different values for these functions to an instance of the
`Simulation` class, you can change how the `Simulation` generates the
phylogenetic tree, simulates sequence evolution, and assigns lengths and
coverages to transcripts.

#### Arguments accepted by custom functions

Importantly, these functions need to be able to accept the kinds of arguments
the `Simulation` will give them. Specifically, all functions need to accept at
least the `py_rand` and `np_rand` arguments, which represent the simulation's
Python PRNG and NumPy PRNG, respectively. (See the [`use_py` and `use_np`
section](#use_py and use_np) for information about two functions that are useful
for making functions that accept these arguments.)

Each of these functions also must accept other positional and keyword arguments 
that are described below.

##### `tree_generator`

| Keyword argument  | Expected type             | Description                                                         |
|-------------------|---------------------------|---------------------------------------------------------------------|
| `taxon_namespace` | `dendropy.TaxonNamespace` | Structure that contains taxa (and their labels) to use in the tree. |
| `num_extant_tips` | `int`                     | Number of extant taxa in the final generated tree.                  |

##### `char_generator`

###### Positional arguments

| Position | Expected type   | Description                              |
|----------|-----------------|------------------------------------------|
| 0        | `int`           | Length of sequence to generate/evolve.   |
| 1        | `dendropy.Tree` | The tree on which to simulate evolution. |

###### Keyword arguments

| Keyword argument | Expected type | Description                                           |
|------------------|---------------|-------------------------------------------------------|
| `root_states`    | `str`         | Sequence with which to begin the simulated evolution. |

Note that the keyword argument above is not used by the base `Simulation` class
but is nevertheless passed to the `char_generator` function by the
`seeded_generate_sequences` function in [`simulation.py`](../simulation.py). For
the `Simulation` class, the value for this keyword argument is always `None`.

##### `seq_len_distribution`

No additional arguments are required.

##### `coverage_distribution`

No additional arguments are required.

#### Registering functions for use in configuration files

You can "register" your custom functions for use in configuration files by
adding keys to the appropriate dictionaries in
`SimulationConfig.load_functions` at runtime. For example, to add a function
called `foo` as a `char_generator`, put the following code before you parse the
configuration file.

```python
SimulationConfig.load_functions["char_generator"]["foo"] = preset(foo)
```

#### Advantages

* Requires adding minimal new code.
* Takes advantage of future updates to the simulation code.
* Makes difference from base `Simulation` clear.

#### Disadvantages

* Limited flexibility.

### Subclassing `Simulation` and `SimulationConfig`

This option is designed to offer a nice balance of flexibility, readability, and
reusability.

You can subclass `Simulation` and `SimulationConfig` and override any of their
methods to customize the behavior of the simulation. You can additionally add
new parts to the simulation by overriding the `simulate` method. See the
[`HybridSimulation`](../hybrid_sim.py) and
[`DivergenceSimulation`](../divergence_sim.py) classes for examples of this
being done.

For using this approach, it is recommended that you read the [Design
overview](#Design overview) above to understand what functions can be overridden
to alter simulation behavior.

#### Advantages

* Somewhat flexible.
* Makes differences from base `Simulation` clear.
* Allows easy reuse of code and updating of the base `Simulation`.

#### Disadvantages

* Composing multiple simulations in arbitrary ways can be difficult.
