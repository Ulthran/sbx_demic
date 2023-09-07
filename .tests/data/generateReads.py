import math
import os
from pathlib import Path
from scipy.stats import cosine


class Genome:
    def __init__(self, fp: Path) -> None:
        self.fp = fp

        with open(self.fp) as f:
            f.readline()  # Annotation line
            self.seq = "".join([l.strip() for l in f.readlines()])

        self.len = len(self.seq)

        self.COMPLEMENT_BASES = {
            "T": "A",
            "C": "G",
            "A": "T",
            "G": "C",
        }

    def sim_reads(self, n: int, a: int, l: int) -> list:
        """
        n: total number of reads to generate
        a: amplitude of sine wave (should end up being PTR value)
        l: length of read
        return: list of generated reads
        """
        idxs = list(cosine.rvs(size=n, scale=1 / a))
        idxs = [int((idx + math.pi) * self.len / (2 * math.pi)) for idx in idxs]
        assert all([idx <= self.len for idx in idxs])
        fwd_reads = [self.get_circular_read(idx, l) for idx in idxs]
        rev_reads = [self.get_rev_read(r) for r in fwd_reads]
        return list(zip(fwd_reads, rev_reads))

    def get_circular_read(self, idx: int, l: int) -> str:
        if idx + l < self.len:
            return self.seq[idx : idx + l]
        else:
            overlap = idx + l - self.len
            print(f"{idx} {l} {self.len}")
            return self.seq[idx:] + self.seq[:overlap]

    def get_rev_read(self, read: str) -> str:
        rc = [self.COMPLEMENT_BASES[x] for x in read]
        rc.reverse()
        return "".join(rc)


class Reads:
    def __init__(self, name: str, reads: list) -> None:
        self.name = name
        self.fwd_reads = [r[0] for r in reads]
        self.rev_reads = [r[1] for r in reads]

    def write_fastqs(self, out_dir: Path):
        self.write_fastq(out_dir / f"{self.name}_R1.fastq", self.fwd_reads)
        self.write_fastq(out_dir / f"{self.name}_R2.fastq", self.rev_reads)

    @staticmethod
    def write_fastq(fp: Path, reads: list):
        with open(fp, "w") as f:
            for i, read in enumerate(reads):
                f.write(f"@{i}\n")
                f.write(f"{read}\n")
                f.write("+\n")
                qual_str = "".join(["~" for r in read])
                f.write(f"{qual_str}\n")


def combine_reads(reads_list: list) -> list:
    return [item for sublist in reads_list for item in sublist]


akk = Genome(Path(os.path.dirname(__file__)) / "reference" / "akk-genome.fasta")
bfrag = Genome(Path(os.path.dirname(__file__)) / "reference" / "Bfragilis.fasta")
ecoli = Genome(Path(os.path.dirname(__file__)) / "reference" / "Ecoli.fasta")

A1 = Reads(
    "A1", combine_reads([akk.sim_reads(10000, 1, 500), bfrag.sim_reads(5000, 2, 500)])
)
A2 = Reads(
    "A2", combine_reads([akk.sim_reads(10000, 5, 500), bfrag.sim_reads(5000, 10, 500)])
)
B1 = Reads(
    "B1",
    combine_reads(
        [
            akk.sim_reads(10000, 1, 500),
            bfrag.sim_reads(5000, 10, 500),
            ecoli.sim_reads(5000, 5, 500),
        ]
    ),
)
B2 = Reads(
    "B2",
    combine_reads(
        [
            akk.sim_reads(10000, 2, 500),
            bfrag.sim_reads(5000, 5, 500),
            ecoli.sim_reads(5000, 5, 500),
        ]
    ),
)
B3 = Reads(
    "B3",
    combine_reads(
        [
            akk.sim_reads(10000, 3, 500),
            bfrag.sim_reads(5000, 1, 500),
            ecoli.sim_reads(5000, 5, 500),
        ]
    ),
)

out_dir = Path(os.path.dirname(__file__)) / "reads"

A1.write_fastqs(out_dir)
A2.write_fastqs(out_dir)
B1.write_fastqs(out_dir)
B2.write_fastqs(out_dir)
B3.write_fastqs(out_dir)
