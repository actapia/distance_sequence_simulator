from divergence_sim import DivergenceSimulationConfig

# This script may be run to perform a simulation using the DivergenceSimulation
# class.
#
# In such simulations, an initial phylogenetic tree is generated, and random
# transcripts are generated and evolved on the phylogenetic tree, producing
# othologs in each extant taxon.
#
# Then, transcripts from some disjoints groups of taxa (of specified size) are
# combined to produce "hybrid" taxa.
#
# These hybrid taxa are then evolved further through as many new, independent
# simulations as there are hybrid taxa. Each of these simulations produces a
# specified number of descendant taxa.

# This is the default configuration to be used when values are not provided by
# the user.
default_config = {
    "count": "match",
    "taxa": "match",
    "log_level": "WARNING",
    "log_transcript_every": 1000,
    "seqid_template": "transcript_{gene}",
    "separate_dirs": False,
    "jobs": 1
}

def main():
    config = DivergenceSimulationConfig(default_config=default_config)
    config.load_config()
    sim = config.make_simulation()
    sim.simulate()

if __name__ == "__main__":
    main()
