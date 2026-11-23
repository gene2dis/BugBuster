#!/bin/sh
# Stand-in for gtdbtk in tests/modules/gtdb_tk_batch.nf.test (bind-mounted over
# the image's /usr/local/bin/gtdbtk by gtdbtk_stand_in.config): no database,
# only the summary files a run leaves. 'classify_wf' writes an ar53 summary
# for every genome of --genome_dir (archaea-only bins: GTDB-Tk writes no bac120
# summary then), or no summary at all when a genome name contains "noclass".
case "$1" in
    --version) echo "gtdbtk: version 2.8.0 Copyright 2017 (stand-in)"; exit 0 ;;
    classify_wf) ;;
    *) echo "gtdbtk stand-in: unsupported command $1" >&2; exit 2 ;;
esac
shift
while [ $# -gt 0 ]; do
    case "$1" in
        --genome_dir) genome_dir="$2"; shift ;;
        --out_dir) out_dir="$2"; shift ;;
    esac
    shift
done
mkdir -p "$out_dir/classify"
echo "stand-in gtdbtk log" > "$out_dir/gtdbtk.log"
ls "$genome_dir" | grep -q noclass && exit 0
{
    printf 'user_genome\tclassification\tclosest_genome_reference\tclosest_genome_reference_radius\tclosest_genome_taxonomy\tclosest_genome_ani\tclosest_genome_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings\n'
    for f in "$genome_dir"/*.fa; do
        g=$(basename "$f" .fa)
        printf '%s\td__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\td__Archaea\ttaxonomic classification defined by topology and ANI\tN/A\tN/A\t90.0\t11\t0.9\tN/A\n' "$g"
    done
} > "$out_dir/classify/gtdbtk.ar53.summary.tsv"
