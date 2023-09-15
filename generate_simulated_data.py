from simulation_config import SimulationConfig

# This script may be run to perform a simulation using the base Simulation
# class.

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
    config = SimulationConfig(default_config=default_config)
    config.load_config()
    sim = config.make_simulation()
    sim.simulate()

if __name__ == "__main__":
    main()
