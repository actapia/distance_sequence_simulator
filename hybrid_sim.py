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
    """A simulation ending with hybridization of taxa.

    This simulation begins with the same process as in the base Simulation
    class. A tree is generated at random, and evolution of randomly generated
    transcript sequeneces is simulated on the tree.

    After the initial simulation, this class selects groups of a specified size
    to hybridize. More specifically, hybridization is simulated by combining all
    the transcripts of the "parent" taxa. Care is taken to name the new
    transcripts in such a way that there are no conflicts, and the origin of
    each transcript in a hybrid may be identified.

    Importantly, the groups selected are disjoint; no taxon is a parent to more
    than one hybrid.

    Attributes:
        count (int):           The number of transcripts to generate per taxon.
        seed (int):            The seed to use for initializing PRNGs.
        tree_generator:        A function used to generate a phylogenetic tree.
        seq_len_distribution:  A function used to select transcript lengths.
        char_generator:        A function that simulates sequence evolution.
        save_tree:             The filename at which to save the generated tree.
        seqid_template (str):  A template string used to assign sequence IDs.
        separate_dirs (bool):  Whether to output files in taxon directories.
        coverage_distribution: A function used to select transcript coverages.
        taxa (int):            The number of extant taxa in the generated tree.
        out_dir:               A path to the initial output directory.
        jobs (int):            The number of parallel jobs to use.
        taxa_prefix (str):     The string that begins each taxon label.
        select (int):          The number of "parent" taxa for each hybrid.
        hybrid_out_dir:        A path to the hybrid simulation output directory.
    """
    def __init__(
            self,
            *args,
            select: Optional[int] = None,
            hybrid_out_dir: Optional[Union[str, Path]] = None,
            **kwargs
    ):
        """Create a HybridSimulation object with the specified settings.
                
        Parameters:
            count (int):           Number of transcripts to generate per taxon.
            seed (int):            Seed for initializing PRNGs.
            tree_generator:        Function that generates phylogenetic trees.
            seq_len_distribution:  Function that generates sequence lengths.
            char_generator:        Function that simulates DNA/RNA evolution.
            save_tree:             Output phylogenetic tree filename.
            seqid_template (str):  Template for naming transcripts.
            separate_dirs (bool):  Whether to put files in separate directories.            
            coverage_distribution: Function that generates transcript coverages.
            taxa (int):            Number of extant taxa in the generated tree.
            out_dir:               Output directory path.
            jobs (int):            Parallel jobs to use.
            taxa_prefix (str):     String that begins each taxon label.
            select (int):          Number of "parent" taxa for each hybrid.
            hybrid_out_dir:        Path to hybrid simulation output directory.
        """
        super().__init__(*args, **kwargs)
        self.select = select
        self.hybrid_out_dir = hybrid_out_dir
    
    def simulate(self) -> tuple[None, dict[str, str]]:
        """Perform the full simulation.

        The simulation begins with the steps from the parent Simulation class.

        Then, (taxa / select) groups of select taxa are randomly chosen to
        hybridize, and their transcripts are combined into new transcript FASTA
        files for the hybrid taxa. These files contain all transcripts from the
        "parent" taxa, but transcripts are renamed to identify their origins.

        Returns:
             A None object. (Since no new phylogenetic tree is generated.)
             A dict mapping taxon labels to paths to transcript FASTA files.
        """
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
    """A configuration manager for creating HybridSimulation objects.

    Attributes:
        default_config (dict):   The default base configuration.
        config (dict):           The current configuration.
        simulation_class (type): The simulation class to construct.
    """
    def __init__(self, *args, simulation_class = HybridSimulation, **kwargs):
        super().__init__(*args, simulation_class=simulation_class, **kwargs)

    def _make_parser(self):
        parser = super()._make_parser()
        parser.add_argument("--hybrid-out-dir", type=Path)
        return parser
    
    def _post_cli(self):
        super()._post_cli()
        self.config.setdefault("select", 2)

        
