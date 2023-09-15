from divergence_sim import DivergenceSimulation, DivergenceSimulationConfig

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
