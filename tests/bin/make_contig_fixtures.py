#!/usr/bin/env python3
"""Generate the committed contig FASTA fixtures in tests/data/contigs/.

Used by tests/modules/bbmap.nf.test (audit #20): a mixed-length set where the
short contig is filtered out, and an all-short set where BBMAP's filter report
must record contigs_after 0 / all_contigs_filtered TRUE (previously the
pipefail-broken 'grep -c || echo 0' counting crashed the task in exactly that
case). Deterministic content - rerunning reproduces the same files.
"""

from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "contigs"


def seq(length: int) -> str:
    s = ("ACGTTGCA" * (length // 8 + 1))[:length]
    return "\n".join(s[i:i + 80] for i in range(0, len(s), 80))


def write_fasta(path: Path, contigs: dict) -> None:
    with open(path, "w") as fh:
        for name, length in contigs.items():
            fh.write(f">{name}\n{seq(length)}\n")
    print(f"wrote {path}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_fasta(OUT_DIR / "mixed_length_contigs.fa",
                {"contig_long1": 1200, "contig_long2": 1500, "contig_short": 300})
    write_fasta(OUT_DIR / "short_only_contigs.fa",
                {"contig_short1": 300, "contig_short2": 250})


if __name__ == "__main__":
    main()
