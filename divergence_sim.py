import re

from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from dendropy import Tree
from joblib import Parallel, delayed
from tqdm import tqdm

from simulation import Simulation, seeded_generate_sequences
from hybrid_sim import HybridSimulation, HybridSimulationConfig
from length_distributions import constant

class RootStateSimulation(Simulation):
    def __init__(
            self,
            *args,
            base_sequences: Optional[Sequence[SeqRecord]] = None,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.base_sequences = base_sequences
        self._seq_id_regex = None

    @property
    def count(self):
        if self._count is None and self.base_sequences is not None:
            self._count = len(self.base_sequences)
        return self._count

    @count.setter
    def count(self, count):
        self._count = count

    @property
    def seq_id_regex(self):
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
    ):
        seeds = self._make_seeds(count)
        return Parallel(n_jobs=jobs, backend="loky", return_as="generator")(
            delayed(
                seeded_generate_sequences
            )(self.char_generator, tree, constant(len(seq.seq)), seed, seq.seq)
            for seed, seq in zip(tqdm(seeds), self.base_sequences)
        )

    def _name_sequence(self, i: int, mat) -> str:
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
    def __init__(
            self,
            *args,
            diverged_taxa: Optional[int] = None,
            divergence_out_dir: Optional[Union[str, Path]] = None,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.diverged_taxa = diverged_taxa
        self.divergence_out_dir = divergence_out_dir

    def _make_child_simulation(self, save_tree, out_dir, prefix, base_sequences):
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
    def __init__(self, *args, simulation_class = DivergenceSimulation, **kwargs):
        super().__init__(*args, simulation_class=simulation_class, **kwargs)
        
    def _make_parser(self):
        parser = super()._make_parser()
        parser.add_argument("--divergence-out-dir", type=Path)
        parser.add_argument("--divergence-out-root", type=Path)
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
