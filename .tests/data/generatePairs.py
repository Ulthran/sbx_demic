import os
from more_itertools import grouper
from pathlib import Path

COMPLEMENT_BASES = {
    "T": "A",
    "C": "G",
    "A": "T",
    "G": "C",
}


def parse_fastq(f):
    for g in grouper(f.readlines(), 4):
        header_str = g[0][1:].strip()
        seq_str = g[1].strip()
        plus_str = g[2].strip()
        quality_str = g[3].strip()

        yield (header_str, seq_str, plus_str, quality_str)


def write_fastq(record, f):
    for i, l in enumerate(record):
        if i == 0:
            f.write(f"@{l}\n")
        else:
            f.write(f"{l}\n")


def write_many_fastq(record_list, f):
    record_list = [
        [f"@{r[0]}\n", f"{r[1]}\n", f"{r[2]}\n", f"{r[3]}\n"] for r in record_list
    ]
    record_list = [item for sublist in record_list for item in sublist]
    f.writelines(record_list)


def get_rev_read(read: str) -> str:
    rc = [COMPLEMENT_BASES[x] for x in read]
    rc.reverse()
    return "".join(rc)


for fn in os.listdir("someReads/"):
    fp = Path("someReads") / fn
    pair_fp = Path("someReads") / fn.replace("_1.fastq", "_2.fastq")
    with open(fp) as f1, open(pair_fp, "w") as f2:
        for h, s, p, q in parse_fastq(f1):
            write_fastq((h, get_rev_read(s), p, q[::-1]), f2)
