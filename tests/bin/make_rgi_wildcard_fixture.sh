#!/bin/bash
#
# Regenerate tests/data/rgi/wildcard_mini/ (used by
# tests/modules/rgi_load_wildcard.nf.test): a synthetic CARD + WildCARD pair
# small enough for a real RGI_LOAD_WILDCARD run.
#   card/card.json                   two protein homolog models in CARD's
#                                    card.json layout (version 1.0.0)
#   card/card_database_v1.0.0.fasta  rgi card_annotation of it, and its
#   card/card_database_v1.0.0_all.fasta  _all twin (what --custom_rgi_card_db
#                                    holds after a first CARD load)
#   card_unannotated/card.json       the same card.json without the
#                                    annotation FASTAs
#   wildcard/                        index-for-model-sequences.txt and the
#                                    eight *_variants.fasta files of a WildCARD
#                                    download (two variants per model in the
#                                    protein homolog files, the others empty),
#                                    no k-mer files
#
# Nothing here is CARD data (CARD's terms exclude commercial use, so the
# public repository carries none): names, accessions, descriptions and
# sequences are made up (seeded random), only the file layouts follow CARD
# 4.0.1 and WildCARD as RGI 6.0.8 reads them. Byte-stable: fixed seed, sorted
# JSON keys, and rgi card_annotation writes the FASTAs in model order.
# Requires docker and python3; no network.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
OUT_DIR="${REPO_DIR}/tests/data/rgi/wildcard_mini"
RGI_IMG=$(grep -o "'quay.io/biocontainers/rgi:[^']*'" "${REPO_DIR}/modules/local/rgi_load_wildcard/main.nf" | tr -d "'")

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}/card" "${OUT_DIR}/card_unannotated" "${OUT_DIR}/wildcard"

python3 -I - "${OUT_DIR}" <<'EOF'
import json, random, sys

out = sys.argv[1]
rng = random.Random(21)
CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a + b + c] = aas[i]
            i += 1
SENSE = [c for c, aa in CODON.items() if aa != "*"]


def gene(n_codons):
    return "ATG" + "".join(rng.choice(SENSE) for _ in range(n_codons)) + "TAA"


def translate(dna):
    return "".join(CODON[dna[i:i + 3]] for i in range(0, len(dna) - 3, 3))


def mutate(dna, n):
    """n synonymous-or-not substitutions inside the coding part, never creating a stop codon."""
    s = list(dna)
    done = 0
    while done < n:
        p = rng.randrange(3, len(s) - 3)
        old = s[p]
        s[p] = rng.choice([b for b in "ACGT" if b != old])
        cod = p - p % 3
        if CODON["".join(s[cod:cod + 3])] == "*":
            s[p] = old
            continue
        done += 1
    return "".join(s)


def category(acc, cvterm, name, cls):
    return {"category_aro_accession": acc, "category_aro_cvterm_id": cvterm, "category_aro_name": name,
            "category_aro_description": f"Synthetic test category {name}.", "category_aro_class_name": cls}


models, index_rows, nuc, prot = {}, [], [], []
specs = [("9001", "testGeneA", "9900001", "test gene family A", 120),
         ("9002", "testGeneB", "9900002", "test gene family B", 150)]
prevalence_id = 0
for k, (model_id, name, aro, family, n_codons) in enumerate(specs):
    dna = gene(n_codons)
    models[model_id] = {
        "model_id": model_id, "model_name": name, "model_type": "protein homolog model",
        "model_type_id": "40292", "model_description": "Synthetic protein homolog model for tests.",
        "model_param": {"blastp_bit_score": {"param_type": "BLASTP bit-score",
                                             "param_description": "Synthetic bit-score cutoff.",
                                             "param_type_id": "40725", "param_value": "150"}},
        "model_sequences": {"sequence": {str(9100 + k): {
            "protein_sequence": {"accession": f"TEST{k + 1:04d}.1", "sequence": translate(dna)},
            "dna_sequence": {"accession": f"TESTN{k + 1:04d}.1", "fmin": "0", "fmax": str(len(dna)),
                             "strand": "+", "sequence": dna, "partial": "0"},
            "NCBI_taxonomy": {"NCBI_taxonomy_cvterm_id": "99001", "NCBI_taxonomy_name": "Testus syntheticus",
                              "NCBI_taxonomy_id": "99001"}}}},
        "ARO_accession": aro, "ARO_id": str(99000 + k), "ARO_name": name, "CARD_short_name": name,
        "ARO_description": f"Synthetic test gene {name}.",
        "ARO_category": {
            str(99100 + k): category(str(9910000 + k), str(99100 + k), family, "AMR Gene Family"),
            "99200": category("9920000", "99200", "testmycin", "Drug Class"),
            "99300": category("9930000", "99300", "test inactivation", "Resistance Mechanism"),
        },
    }
    for v in range(2):
        prevalence_id += 1
        var = mutate(dna, 3 + v)
        head = (f"Prevalence_Sequence_ID:{prevalence_id}|ARO_Name:{name}|ARO:{aro}|"
                f"Detection_Model:Protein Homolog Model|CARD_Short_Name:{name}")
        nuc.append(f">{head}\n{var}\n")
        prot.append(f">{head}\n{translate(var)}\n")
        index_rows.append([str(prevalence_id), model_id, name, f"ARO:{aro}", "protein homolog model",
                           "Testus syntheticus", f"TESTC{prevalence_id:04d}", "ncbi_contig", "Strict", "99.0",
                           "300.0", family, "test inactivation", "testmycin", name])

card = dict(models)
card["_version"] = "1.0.0"
card["_comment"] = {"description": "Synthetic CARD-layout test file (BugBuster tests); no CARD data.",
                    "access": "public"}
card["_timestamp"] = "2026-11-22T00:00:00+00:00"
for d in ("card", "card_unannotated"):
    with open(f"{out}/{d}/card.json", "w") as f:
        json.dump(card, f, sort_keys=True)

header = ["prevalence_sequence_id", "model_id", "aro_term", "aro_accession", "detection_model", "species_name",
          "ncbi_accession", "data_type", "rgi_criteria", "percent_identity", "bitscore", "amr_gene_family",
          "resistance_mechanism", "drug_class", "card_short_name"]
with open(f"{out}/wildcard/index-for-model-sequences.txt", "w") as f:
    f.write("\t".join(header) + "\n")
    f.writelines("\t".join(r) + "\n" for r in index_rows)
for kind in ("protein_homolog", "protein_overexpression", "protein_variant", "rRNA_gene_variant"):
    with open(f"{out}/wildcard/nucleotide_fasta_{kind}_model_variants.fasta", "w") as f:
        f.writelines(nuc if kind == "protein_homolog" else [])
    if kind != "rRNA_gene_variant":
        with open(f"{out}/wildcard/protein_fasta_{kind}_model_variants.fasta", "w") as f:
            f.writelines(prot if kind == "protein_homolog" else [])
EOF

docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${OUT_DIR}/card:/card" -w /card "${RGI_IMG}" \
    bash -c 'rgi card_annotation -i card.json > /tmp/card_annotation.log 2>&1 || { cat /tmp/card_annotation.log >&2; exit 1; }'

cat > "${OUT_DIR}/README.md" <<'EOF'
Synthetic CARD + WildCARD fixture for tests/modules/rgi_load_wildcard.nf.test,
written by tests/bin/make_rgi_wildcard_fixture.sh (seeded, byte-stable). It
follows the layouts of CARD's card.json and a WildCARD variants download but
holds no CARD data: every name, accession, description and sequence is made
up.
EOF
ls -la "${OUT_DIR}"/*
