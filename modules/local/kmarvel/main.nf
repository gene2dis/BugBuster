/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KMARVEL Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    K-MARVEL read-level ARG detection: translated k-mer search of the reads
    against the CARD protein models, gene-centric protein assembly, chimera
    removal and SNP-aware variant-model filtering, then a blastp check of
    every assembled ARG (https://bitbucket.org/amr-avengers/k-marvel).

    Container: K-MARVEL has no bioconda/biocontainer package. Every engine
    uses the pipeline-built image (containers/kmarvel/Dockerfile, published
    to GHCR by .github/workflows/kmarvel-image.yml); pin it by digest after
    each rebuild.

    Upstream behaviours handled here (commit 7108d2d):
      - it reads gzip input with flate2's single-member GzDecoder, so a `cat`
        of several .gz files (and fastp's multi-member output) is silently
        truncated after the first member: the reads are decompressed into
        one plain FASTQ instead
      - output and temporary file names are cut at the first '.', so it runs
        on fixed dot-free names and the results are renamed to the sample
      - a failed blastp exits 0, and zero ARGs give a header-less TSV:
        bin/kmarvel_check_output.sh handles both
      - --reference is also the blastp -db, so the FASTA is copied into the
        task and indexed there (the shared database is never written to)
      - -k/-e are ignored (fixed by --data-type) and are not exposed

    Input:
        tuple val(meta), path(reads)
        reads: host-removed reads [R1, R2, optional singletons]; all are
               concatenated, as upstream asks for paired-end data
        kmarvel_db: CARD_protein_all.fasta, metadata.tsv, snps.txt
                    (+ DB_VERSION when built by FORMAT_KMARVEL_DB)

    Output:
        tsv: tuple val(meta), path("*_kmarvel.tsv")
        json: path("*_kmarvel.json")
        log: path("*_kmarvel.log")
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process KMARVEL {
    tag "${meta.id}"
    // 1.0.0-0fd49c9 (containers/kmarvel); re-pin after a rebuild
    container 'ghcr.io/gene2dis/bugbuster-kmarvel@sha256:0fd49c9a1665d63541913f090ebc8cbaa75d4e60753addaf6ff34f5e2747b5a5'

    label 'process_medium'

    input:
        tuple val(meta), path(reads)
        path(kmarvel_db)

    output:
        tuple val(meta), path("*_kmarvel.tsv"), emit: tsv
        path("*_kmarvel.json"), emit: json
        path("*_kmarvel.log"), emit: log
        path("versions.yml"), emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        """
        for f in CARD_protein_all.fasta metadata.tsv snps.txt; do
            if [ ! -s ${kmarvel_db}/\$f ]; then
                echo "ERROR: ${kmarvel_db} has no \$f (expected a FORMAT_KMARVEL_DB / --custom_kmarvel_db layout)" >&2
                exit 1
            fi
        done

        # Own copy of the reference with its BLAST index; the small metadata
        # files are only read
        cp -L ${kmarvel_db}/CARD_protein_all.fasta CARD_protein_all.fasta
        makeblastdb -in CARD_protein_all.fasta -dbtype prot -blastdb_version 5 \\
            -title CARD_protein_all -out CARD_protein_all.fasta > makeblastdb.log
        ln -s ${kmarvel_db}/metadata.tsv metadata.tsv
        ln -s ${kmarvel_db}/snps.txt snps.txt

        zcat ${reads} > reads.fastq

        kmarvel \\
            --filename reads.fastq \\
            --output kmarvel.json \\
            --data-type sr \\
            --reference CARD_protein_all.fasta \\
            --metadata metadata.tsv \\
            --snpdata snps.txt \\
            --threads ${task.cpus} \\
            ${args} \\
            2>&1 | tee kmarvel.log

        rm -f reads.fastq

        kmarvel_check_output.sh kmarvel.log kmarvel.tsv

        mv kmarvel.json ${prefix}_kmarvel.json
        mv kmarvel.tsv ${prefix}_kmarvel.tsv
        mv kmarvel.log ${prefix}_kmarvel.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kmarvel: \$(kmarvel --version 2>&1 | sed 's/^kmarvel //')
            blast: \$(blastp -version 2>&1 | head -n 1 | sed 's/^blastp: //; s/+\$//')
            card: "\$(cat ${kmarvel_db}/DB_VERSION 2>/dev/null || echo custom)"
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"
        """
        printf 'aro_full_name\\taro_id\\tgene_family\\tdrug_class\\tresistance_mechanism\\tmodel\\tidentity\\tassembly_coverage\\tsupporting_reads\\n' > ${prefix}_kmarvel.tsv
        echo '{}' > ${prefix}_kmarvel.json
        touch ${prefix}_kmarvel.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kmarvel: 1.0.0
            blast: 2.17.0
            card: stub
        END_VERSIONS
        """
}
