import contextlib
import logging
import random

import numpy as np
import dendropy.utility

from pathlib import Path
from typing import Callable, Optional, Union, TextIO, Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dendropy import Tree, Taxon
from joblib import Parallel, delayed

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda x: x

def seeded_generate_sequences(generate_chars, tree, seq_len_dist, seed):
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
        np_rand=np_rand
    )
    return list(res.values())

class Simulation:
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
            jobs: Optional[int] = None
    ):
        self.count = count
        self.seed = seed
        self.tree_generator = tree_generator
        self.seq_len_distribution = seq_len_distribution
        self.char_generator = char_generator
        self.save_tree = save_tree
        self.seqid_template = seqid_template
        self.separate_dirs = separate_dirs
        self.coverage_distribution = coverage_distribution

    def rand_init(self):
        self.rand = random.Random()
        # Set random seed if provided.
        try:
            self.rand.seed(self.seed)
            self.np_rand = np.random.RandomState(seed=self.seed+1)
        except AttributeError:
            self.np_rand = np.random.RandomState()

    def _generate_tree(self):
        logging.debug("Generating tree.")
        tree = self.tree_generator(py_rand=self.rand, np_rand=self.np_rand)
        return tree

    def _generate_coverage(self, i, mat):
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
    ):
        seeds = self.np_rand.randint(np.iinfo(np.int32).max, size=count)
        return Parallel(n_jobs=jobs, backend="loky", return_generator=True)(
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

    def simulate(self) -> tuple[Tree, dict[str, str]]:
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
                cov = self._generate_coverage(i, mat)
                # TODO: Maybe allow more generic way of specifying sequence
                # attributes to generate.
                seq_id = self.seqid_template.format(
                    cov=cov,
                    gene=i,
                    # At some point, we may allow multiple isotigs.
                    iso=1
                )
                for f, seq in zip(out_files.values(), mat):
                    SeqIO.write(
                        SeqRecord(Seq(str(seq)), id=seq_id, description=""),
                        f,
                        format="fasta"
                    )
            out_paths = {k: v.name for (k, v) in out_files.items()}
        return tree, out_paths
