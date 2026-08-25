#!/usr/bin/env python3
"""Generate the committed bin-directory fixtures used by
tests/modules/prodigal_bins.nf.test (regression tests for audit item #22).

Creates, under tests/data/bins/:
  - marker_only_bins/SKIPPED.txt  a MetaWRAP-style marker dir with no *.fa
    (what METAWRAP emits when there is nothing to refine)
  - real_bins/bin.1.fa            one 25 kb deterministic pseudo-random contig
    (prodigal -p single refuses sequences under 20 kb)

Deterministic: re-running reproduces byte-identical fixtures.
"""

import random
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "bins"


def main():
    marker_dir = OUT_DIR / "marker_only_bins"
    marker_dir.mkdir(parents=True, exist_ok=True)
    (marker_dir / "SKIPPED.txt").write_text("No bins available for refinement\n")

    real_dir = OUT_DIR / "real_bins"
    real_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)
    seq = "".join(rng.choice("ACGT") for _ in range(25000))
    lines = [seq[i:i + 80] for i in range(0, len(seq), 80)]
    (real_dir / "bin.1.fa").write_text(">bin1_contig_1\n" + "\n".join(lines) + "\n")

    print(f"Fixtures written under {OUT_DIR}")


if __name__ == "__main__":
    main()
