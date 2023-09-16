import contextlib
import logging
import random

import numpy as np
import numpy.typing
import dendropy.utility

from pathlib import Path
from typing import Callable, Optional, Union, TextIO, Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dendropy import Tree, Taxon, TaxonNamespace, DnaCharacterMatrix
from joblib import Parallel, delayed

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda x: x

def seeded_generate_sequences(
        generate_chars: Callable[..., DnaCharacterMatrix],
        tree: Tree,
        seq_len_dist: Callable[..., int],
        seed: int,
        base_sequence: Optional[str] = None
) -> list[str]:
    """Simulate evolution of sequences on the given tree using given RNG seed.

    This function initializes a NumPy and Python pseudo-random number generator
    with the specified seed and then uses the provided seq_len_dist and
    generate_chars functions to simulate random evolution of a random-length
    sequence along a provided tree. (Actually, the NumPy RNG is initialized with
    a seed one greater than the provided seed.)

    Optionally, this function may be provided with a "base sequence" (root
    state) from which to start the evolution. If no such base sequence is
    provided, then a sequence will be randomly generated.

    Parameters:
        generate_chars:      A function to simulate DNA/RNA evolution on a tree.
        tree:                The tree on which to perfrom the simulation.
        seq_len_dist:        A function that returns a random sequence length.
        seed (int):          A seed used to initialize the PRNGs.
        base_sequence (str): The root state of the DNA/RNA sequence.

    Returns:
        A list of sequences for each extant taxon in the tree, in order.
    """
    seed = int(seed)
    rand = random.Random()
    rand.seed(seed)
    # By design, the line below causes this code to fail if the DendroPy version
    # is too old to support reproducible pseudo-random number generation
    # correctly.
    dendropy.utility.GLOBAL_RNG = None
    np_rand = np.random.RandomState(seed=seed+1)
    res = generate_chars(
        seq_len_dist(py_rand=rand, np_rand=np_rand),
        tree,
        py_rand=rand,
        np_rand=np_rand,
        root_states=base_sequence
    )
    return list(res.values())

class Simulation:
    """A simulation in which random sequences are evolved on a generated tree.

    Broadly, the simulation involves two parts.

    First, a random phylogenetic tree is generated using the configured
    tree_generator function. The number of extant taxa in the generated tree is
    configurable.

    Second, a number of random transcripts of random length are generated, and
    their evolution along the tree is simulated. This results in orthologs of
    the transcript for each extant taxon.

    The outputs of the simulation may be configured. By default, the Simulation
    outputs to a specified out_dir a FASTA file containing the transcripts
    generated for each taxon.

    The parts of the simulation are performed in separate functions to
    facilitate overriding their behavior in descendant classes of this
    Simulation class.

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
        out_dir:               A path to the output directory.
        jobs (int):            The number of parallel jobs to use.
        taxa_prefix (str):     The string that begins each taxon label.        
    """
    def __init__(
            self,
            count: Optional[int] = None,
            seed: Optional[int] = None,
            
            # TODO: Specify these types more completely.
            tree_generator: Optional[Callable[..., Tree]] = None,
            seq_len_distribution: Optional[Callable[..., int]] = None,
            char_generator: Optional[Callable] = None,
            
            save_tree: Optional[Union[str, Path]] = None,
            seqid_template: Optional[str] = None,
            separate_dirs: Optional[bool] = None,
            coverage_distribution: Optional[bool] = None,

            taxa: Optional[int] = None,
            out_dir: Optional[Union[str, Path]] = None,
            jobs: Optional[int] = None,

            taxa_prefix: Optional[str] = None
    ):
        """Create a Simulation object with the specified settings.

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
        """
        self.count = count
        self.seed = seed
        self.tree_generator = tree_generator
        self.seq_len_distribution = seq_len_distribution
        self.char_generator = char_generator
        self.save_tree = save_tree
        self.seqid_template = seqid_template
        self.separate_dirs = separate_dirs
        self.coverage_distribution = coverage_distribution
        self.taxa = taxa
        self.out_dir = out_dir
        self.jobs = jobs
        self.taxa_prefix = taxa_prefix

    def _make_seeds(self, count: int) -> np.typing.NDArray[np.int32]:
        """Generate a specified number of integers to be used as PRNG seeds.

        Such seeds are necessary for "child" PRNGs used, e.g., in parallel
        worker processes.

        Parameters:
            count (int): The number of seeds to generate

        Returns:
            A NumPy array containing count integers to be used as PRNG seeds.
        """
        return self.np_rand.randint(np.iinfo(np.int32).max, size=count)

    def _make_seed(self) -> np.int32:
        """Generate an integer to be used as a PRNG seed.

        See the documentation for Simulation._make_seeds for an explanation of
        why this is useful.

        Returns:
            An integer to be used as a PRNG seed.
        """
        return self._make_seeds(1)[0]

    def _make_taxon_namespace(self) -> TaxonNamespace:
        """Create a taxon namespace for this simulation.

        This taxon namespace will be used for assigning taxa to the extant
        tips of the generated phylogenetic tree.

        In this Simulation class, the taxa are simply named by appending
        successive integers to the configured taxa prefix, starting from 0.
        The total number of taxa in the returned namespace is equal to the
        number of extant taxa to be generated in the tree.

        Returns:
            A TaxonNamespace object for the taxa in the generated tree.
        """        
        return TaxonNamespace(
            [f"{self.taxa_prefix}{i}" for i in range(self.taxa)]
        )

    def rand_init(self):
        """Make/initialize the simulation's pseudo-random number generators."""
        self.rand = random.Random()
        # Set random seed if provided.
        try:
            self.rand.seed(self.seed)
            self.np_rand = np.random.RandomState(seed=self.seed+1)
        except AttributeError:
            self.np_rand = np.random.RandomState()

    def _generate_tree(self) -> Tree:
        """Create a phylogenetic tree using the configured tree generator."""
        logging.debug("Generating tree.")
        ns = self._make_taxon_namespace()
        try:
            tree = self.tree_generator(
                py_rand=self.rand,
                np_rand=self.np_rand,
                taxon_namespace=ns,
                num_extant_tips=self.taxa
            )
        except TypeError:
            # This except clause is present in case the num_extant_tips argument
            # was already specified.
            tree = self.tree_generator(
                py_rand=self.rand,
                np_rand=self.np_rand,
                taxon_namespace=ns,
            )
        return tree

    def _generate_coverage(self, i: int, mat: DnaCharacterMatrix) -> float:
        """Assign a coverage value for the given transcript.

        This function uses the configured coverage distribution to generate a
        (possibly random) coverage value for a transcript.

        This function accepts the index and DnaCharacterMatrix for the
        transcript to which the coverage is being assigned. In this class, these
        values are currently unused, but they are present in case a subclass
        needs to base the coverage value on the index or sequence of the
        transcript. These parameters may also be used for some coverage
        distributions in a future version of this code.

        Parameters:
            i (int): The index of the transcript.
            mat:     The generated DnaCharacterMatrix for the transcript.

        Returns:
            A coverage value for the specified transcript.
        """
        try:
            return self.coverage_distribution(
                py_rand=self.rand,
                np_rand=self.np_rand
            )
        except AttributeError:
            return float(self.count - i)

    def _generate_sequences(
            self,
            tree: Tree,
            count: int = 1,
            jobs: int = 1
    ) -> Iterable[list[str]]:
        """Generate transcript sequences for each extant taxon on the tree.

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
            )(self.char_generator, tree, self.seq_len_distribution, seed)
            for seed in tqdm(seeds)
        )

    def _open_out_files(
            self,
            stack: contextlib.ExitStack,
            taxa : Iterable[Taxon],
    ) -> dict[str, TextIO]:
        """Open output FASTA files on the stack for all taxa.

        Parameters:
            stack: The ExitStack on which to open the files.
            taxa:  The taxa for which to open output files.

        Returns:
            A dictionary mapping taxa labels to file-like objects.
        """
        logging.debug("Opening output files for writing.")
        out_files = {}
        for t in taxa:
            if self.separate_dirs:
                taxon_dir = self.out_dir / Path(t.label)
                taxon_dir.mkdir(exist_ok=True)
                transcripts_file = taxon_dir / "transcripts.fasta"
                f = stack.enter_context(open(transcripts_file, "w"))
            else:
                transcripts_file = self.out_dir / Path(
                    "{}.fasta".format(t.label)
                )
                f = stack.enter_context(
                    open(
                        transcripts_file,
                        "w"
                    )
                )
            out_files[t.label] = f
        return out_files

    def _name_sequence(self, i: int, mat) -> str:
        """Assign a name to a transcript.

        This function uses the configured name template to generate a name.

        The function accepts the index and DnaCharacterMatrix for the transcript
        to which the name is being assigned. In this class, these values are
        currently unused, but they are present for potential use by subclasses.

        Parameters:
            i (int): The index of the transcript.
            mat:     The generated DnaCharacterMatrix for the transcript.

        Returns:
            A name for the specified transcript.        
        """
        cov = self._generate_coverage(i, mat)
        # TODO: Maybe allow more generic way of specifying sequence
        # attributes to generate.
        return self.seqid_template.format(
            cov=cov,
            gene=i,
            # At some point, we may allow multiple isotigs.
            iso=1
        )

    def simulate(self) -> tuple[Tree, dict[str, str]]:
        """Perform the full simulation.

        This function generates a tree, generate and simulates evolution of
        sequences on that tree, and outputs the results to separate FASTA files
        for each taxon.

        For use in "client code" and subclasses, this function returns the
        generated Tree and a dictionary mapping taxon labels to paths to their
        corresponding FASTA files.

        Returns:
            The generated tree.
            A dict mapping taxon labels to paths to transcript FASTA files.
        """
        self.rand_init()
        self.out_dir.mkdir(exist_ok=True)
        tree = self._generate_tree()
        logging.debug(tree.as_ascii_plot())
        if self.save_tree:
            with open(
                    Path(self.out_dir)/self.save_tree,
                    "w"
            ) as out_file:
                out_file.write(tree.as_string(schema="newick"))
        with contextlib.ExitStack() as stack:
            out_files = self._open_out_files(stack, tree.taxon_namespace)
            for i, mat in enumerate(
                    self._generate_sequences(
                        tree,
                        self.count,
                        jobs=self.jobs
                    )
            ):
                # if i%self.config["log_transcript_every"] == 0:
                #     logging.info(f"Writing transcript {i}")
                seq_id = self._name_sequence(i, mat)
                for f, seq in zip(out_files.values(), mat):
                    SeqIO.write(
                        SeqRecord(Seq(str(seq)), id=seq_id, description=""),
                        f,
                        format="fasta"
                    )
            out_paths = {k: v.name for (k, v) in out_files.items()}
        return tree, out_paths
