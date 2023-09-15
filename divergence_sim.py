import re

from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Union, Iterable

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from dendropy import Tree
from joblib import Parallel, delayed
from tqdm import tqdm

from simulation import Simulation, seeded_generate_sequences
from hybrid_sim import HybridSimulation, HybridSimulationConfig
from length_distributions import constant

class RootStateSimulation(Simulation):
    """A simulation that evolves existing transcripts.

    This simulation operates similarly to the base Simulation class, but instead
    of generating new sequences for each transcript, the simulation evolves
    existing sequences that are provided as the "root state" or "base sequences"
    for the sequence evolution simulation.

    Attributes:
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
        base_sequences:        The sequences to use as root states.
    """
    def __init__(
            self,
            *args,
            base_sequences: Optional[Sequence[SeqRecord]] = None,
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
            base_sequences:        Sequences to use as root states.
        """
        super().__init__(*args, **kwargs)
        self.base_sequences = base_sequences
        self._seq_id_regex = None

    @property
    def count(self) -> int:
        """Returns the number of sequences to generate.

        For this class, the count should be the same as the number of
        base_sequences. If the count was not set manually, it will be calculated
        automatically from the base_sequences attribute.
        """
        if self._count is None and self.base_sequences is not None:
            self._count = len(self.base_sequences)
        return self._count

    @count.setter
    def count(self, count: int):
        self._count = count

    @property
    def seq_id_regex(self) -> re.Pattern:
        """Returns a regular expression for parsing the base sequence IDs.

        The regular expression is built using the seqid_template attribute. Once
        the regular expression has been compiled, it is cached.
        """
        if self._seq_id_regex is None:
            self._seq_id_regex = re.compile(
                self.seqid_template.format(
                    **{p: f"(?P<{p}>.*)" for p in ["cov", "gene", "iso"]}
                )
            )
        return self._seq_id_regex

    def _generate_sequences(
            self,
            tree: Tree,
            count: int = 1,
            jobs: int = 1
    ) -> Iterable[list[str]]:
        """Generate transcript sequences for each extant taxon on the tree.

        Unlike _generate_sequences in the parent Simulation class, this function
        does not generate sequences from scratch. Instead, it uses the
        base_sequences provided to this RootStateSimulation object.

        Parameters:
            tree:        The tree on which to simulate sequence evolution.
            count (int): The number of sequences to generate per taxon.
            jobs (int):  The number of parallel jobs to run.

        Returns:
            An iterable of count lists of transcripts for each taxon.
        """

        seeds = self._make_seeds(count)
        return Parallel(n_jobs=jobs, backend="loky", return_as="generator")(
            delayed(
                seeded_generate_sequences
            )(self.char_generator, tree, constant(len(seq.seq)), seed, seq.seq)
            for seed, seq in zip(tqdm(seeds), self.base_sequences)
        )

    def _name_sequence(self, i: int, mat) -> str:
        """Assign a name to a transcript.

        This function generates a name for a transcript based on the
        corresponding base transcript's name. The only modification made is that
        the coverage is altered.

        Parameters:
            i (int): The index of the transcript.
            mat:     The generated DnaCharacterMatrix for the transcript.

        Returns:
            A name for the specified transcript.        
        """
        try:
            cov_span = self.seq_id_regex.search(
                self.base_sequences[i].id
            ).span("cov")
        except AttributeError:
            from IPython import embed
            embed()
        cov = self._generate_coverage(i, mat)
        return "{}{}{}".format(
            self.base_sequences[i].id[:cov_span[0]],
            cov,
            self.base_sequences[i].id[cov_span[1]:]
        )
            

class DivergenceSimulation(HybridSimulation):
    """A simulation ending with hybridization and divergences of taxa.

    This simulation uses the output from a HybridSimulation as input for a next
    step of the simulation in which the sequences from the hybrids are evolved
    along new randomly generataed trees to create transcripts for descendent
    taxa.

    Currently, each new tree must have the same number of extant tips. The new
    trees may represent different total amounts of time; the distances from
    their roots to their tips may differ.

    Attributes:
        count (int):           The number of transcripts to generate per taxon.
        seed (int):            The seed to use for initializing PRNGs.
        tree_generator:        A function used to generate phylogenetic trees.
        seq_len_distribution:  A function used to select transcript lengths.
        char_generator:        A function that simulates sequence evolution.
        save_tree:             The filename at which to save generated trees.
        seqid_template (str):  A template string used to assign sequence IDs.
        separate_dirs (bool):  Whether to output files in taxon directories.
        coverage_distribution: A function used to select transcript coverages.
        taxa (int):            The number of extant taxa in the original tree.
        out_dir:               A path to the initial output directory.
        jobs (int):            The number of parallel jobs to use.
        taxa_prefix (str):     The string that begins each taxon label.
        select (int):          The number of "parent" taxa for each hybrid.
        hybrid_out_dir:        A path to the hybrid simulation output directory.
        diverged_taxa (int):   The number of extant taxa in the new trees.
        divergence_out_dir:    A path to the divergence simulation output dir.
    """
    def __init__(
            self,
            *args,
            diverged_taxa: Optional[int] = None,
            divergence_out_dir: Optional[Union[str, Path]] = None,
            **kwargs
    ):
        """Create a DivergenceSimulation object with the specified settings.
                
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
            taxa (int):            Number of extant taxa in original tree.
            out_dir:               Output directory path.
            jobs (int):            Parallel jobs to use.
            taxa_prefix (str):     String that begins each taxon label.
            select (int):          Number of "parent" taxa for each hybrid.
            hybrid_out_dir:        Path to hybrid simulation output directory.
            diverged_taxa (int):   Number of extant taxa in new trees.
            divergence_out_dir:    Path to the divergence simulation output dir.
        """
 
        super().__init__(*args, **kwargs)
        self.diverged_taxa = diverged_taxa
        self.divergence_out_dir = divergence_out_dir

    def _make_child_simulation(
            self,
            save_tree: str,
            out_dir: Union[str, Path],
            prefix: str,
            base_sequences: Sequence[SeqRecord]
    ) -> RootStateSimulation:
        """Create a RootStateSimulation to simulate evolution of hybrids.

        This function creates a "child" simulation that may be used by this
        simulation to simulate further evolution of some taxon.

        Parameters:
            save_tree (str): Output phylogenetic tree filename.
            out_dir:         Output directory for child simulation.
            prefix:          String that begins taxa labels in child simulation.
            base_sequences:  Sequences to evolve/use as root states.

        Returns:
            A RootStateSimulation configured as specified.
        """
        return RootStateSimulation(
            seed = int(self._make_seed()),
            tree_generator = self.tree_generator,
            char_generator = self.char_generator,
            save_tree = save_tree,
            seqid_template = self.seqid_template,
            separate_dirs = self.separate_dirs,
            coverage_distribution = self.coverage_distribution,
            taxa = self.diverged_taxa,
            out_dir = out_dir,
            jobs = self.jobs,
            taxa_prefix = prefix,
            base_sequences = base_sequences
        )

    def simulate(self):
        """Perform the full simulation.

        The simulation begins with the steps from the parent HybridSimulation
        class.

        Evolution is then simulated for each of the resulting hybrid taxa. For
        each hybrid taxon, this generates four new descendent taxa.

        Returns:
            The generated tree.
            A dict mapping taxon labels to paths to transcript FASTA files.
        """
 
        _, paths = super().simulate()
        self.divergence_out_dir.mkdir(exist_ok=True)
        trees = []
        out_paths = {}
        for taxon_name, base_sequences_path in paths.items():
            if self.separate_dirs:
                out_dir = self.divergence_out_dir / taxon_name
                out_dir.mkdir(exist_ok=True)
                save_tree = self.save_tree
            else:
                out_dir = self.divergence_out_dir
                try:
                    save_tree = f"{taxon_name}.tree"
                except AttributeError:
                    save_tree = None
            child_sim = self._make_child_simulation(
                save_tree,
                out_dir,
                taxon_name + "_T",
                list(SeqIO.parse(base_sequences_path, "fasta"))
            )
            new_tree, new_out_paths = child_sim.simulate()
            trees.append(new_tree)
            out_paths.update(new_out_paths)
        return trees, out_paths

class DivergenceSimulationConfig(HybridSimulationConfig):
    """A configuration manager for creating DivergenceSimulation objects.

    Attributes:
        default_config (dict):   The default base configuration.
        config (dict):           The current configuration.
        simulation_class (type): The simulation class to construct.
    """
    def __init__(self, *args, simulation_class = DivergenceSimulation, **kwargs):
        super().__init__(*args, simulation_class=simulation_class, **kwargs)
        
    def _make_parser(self):
        parser = super()._make_parser()
        parser.add_argument("--divergence-out-dir", type=Path)
        # If divergence_out_root is specified (whether via CLI or the YAML
        # config), then the out_dir, hybrid_out_dir, and divergence_out_dir
        # directories can be created automatically as child directories of
        # divergence_out_root.
        parser.add_argument(
            "--divergence-out-root",
            type=Path,
            help="root directory at which to create output directories"
        )
        return parser

    def _post_cli(self):
        super()._post_cli()
        try:
            self.config["divergence_out_root"].mkdir(exist_ok=True)
        except (KeyError, TypeError):
            pass
        self.config.setdefault(
            "out_dir",
            self.config["divergence_out_root"] / "p1"
        )
        self.config.setdefault(
            "hybrid_out_dir",
            self.config["divergence_out_root"] / "p2"
        )
        self.config.setdefault(
            "divergence_out_dir",
            self.config["divergence_out_root"] / "p3"
        )
