from pathlib import Path
from typing import Optional, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from more_itertools import chunked

from simulation import Simulation
from simulation_config import SimulationConfig

# TODO: Come up with a more modular way of building simulation pipelines.
# For example, HybridSimulation should be something we can apply to the output
# of ANY kind of simulation, not just the base Simulation class.

class HybridSimulation(Simulation):
    def __init__(
            self,
            *args,
            select: Optional[int] = None,
            hybrid_out_dir: Optional[Union[str, Path]] = None,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.select = select
        self.hybrid_out_dir = hybrid_out_dir
    
    def simulate(self):
        if self.taxa % self.select:
            raise ValueError(
                "Number of inital extant taxa must be 0 mod select."
            )
        tree, paths = super().simulate()
        # Further simulation ...
        # Determine which taxa should be hybridized.
        taxa = [t.label for t in tree.taxon_namespace]
        self.rand.shuffle(taxa)
        self.hybrid_out_dir.mkdir(exist_ok=True)
        out_paths = {}
        for i, hybridized in enumerate(chunked(taxa, self.select)):
            taxon_label = "+".join(hybridized)
            if self.separate_dirs:
                taxon_dir = self.hybrid_out_dir / Path(taxon_label)
                taxon_dir.mkdir(exist_ok=True)
                transcripts_path = taxon_dir / "transcripts.fasta"
            else:
                transcripts_path = self.hybrid_out_dir / f"{taxon_label}.fasta"
            out_paths[taxon_label] = transcripts_path
            with open(transcripts_path, "w") as out_file:
                for transcripts in zip(
                        *(
                            SeqIO.parse(paths[taxon], "fasta")
                            for taxon in hybridized
                        )
                ):
                    for taxon, transcript in zip(hybridized, transcripts):
                        SeqIO.write(
                            SeqRecord(
                                transcript.seq,
                                id="{}_{}".format(transcript.id, taxon),
                                description=""
                            ),
                            out_file,
                            format="fasta"
                        )
        return None, out_paths
    
class HybridSimulationConfig(SimulationConfig):
    def __init__(self, *args, simulation_class = HybridSimulation, **kwargs):
        super().__init__(*args, simulation_class=simulation_class, **kwargs)

    def _make_parser(self):
        parser = super()._make_parser()
        parser.add_argument("--hybrid-out-dir", type=Path)
        return parser
    
    def _post_cli(self):
        super()._post_cli()
        self.config.setdefault("select", 2)

        
