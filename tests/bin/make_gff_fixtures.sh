#!/bin/bash
#
# Generate the committed gene-GFF fixtures used by
# tests/modules/featurecounts_genes.nf.test and tests/bin/test_gff2saf.sh
# (functional annotation T3: featureCounts over Pyrodigal gene coordinates).
#
# Creates, under tests/data/gff/:
#   - test_contig1_genes.gff.gz  Pyrodigal-format GFF3 for contig_1 (the 1 kb
#     contig in tests/data/bam/test_bin.1_all_reads.bam): two CDS rows with
#     Prodigal-style ID attributes (ID=1_1 at 1-300 overlapping the fixture
#     reads, ID=1_2 at 601-900 with no reads)
#   - test_empty_genes.gff.gz    header-only GFF (empty-contig sample shape:
#     no genes predicted) -> gff2saf.sh emits a header-only SAF and
#     FEATURECOUNTS_GENES takes its short-circuit path
#
# The ID attributes deliberately use Prodigal's <seqnum>_<genenum> form so the
# tests exercise the GeneID rewrite to <contig>_<genenum> (design doc 4.4).
#
# Requirements: gzip only.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUT_DIR="${SCRIPT_DIR}/../data/gff"
mkdir -p "${OUT_DIR}"
OUT_DIR="$( cd "${OUT_DIR}" && pwd )"

cat > "${OUT_DIR}/test_contig1_genes.gff" <<'EOF'
##gff-version  3
# Sequence Data: seqnum=1;seqlen=1000;seqhdr="contig_1"
# Model Data: version=pyrodigal.v3.6.3;run_type=Metagenomic;model="5|Mycoplasma_bovis_PG45|B|29.3|11|24"
contig_1	pyrodigal_v3.6.3	CDS	1	300	31.6	+	0	ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.500;conf=99.99;score=31.60;cscore=30.00;sscore=1.60;rscore=0.00;uscore=0.00;tscore=1.60;
contig_1	pyrodigal_v3.6.3	CDS	601	900	25.2	-	0	ID=1_2;partial=00;start_type=ATG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;gc_cont=0.480;conf=99.90;score=25.20;cscore=22.00;sscore=3.20;rscore=2.00;uscore=0.40;tscore=0.80;
EOF

cat > "${OUT_DIR}/test_empty_genes.gff" <<'EOF'
##gff-version  3
# Sequence Data: seqnum=1;seqlen=0;seqhdr="empty"
EOF

# -n keeps the gzip output byte-stable (no mtime/name in the header)
gzip -nf "${OUT_DIR}/test_contig1_genes.gff" "${OUT_DIR}/test_empty_genes.gff"

echo "Fixtures written under ${OUT_DIR}"
