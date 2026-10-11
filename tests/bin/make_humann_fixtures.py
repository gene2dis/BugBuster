#!/usr/bin/env python3
"""Build the deterministic HUMAnN read-branch fixture (design doc T8c,
Section 4.6.1). Called twice by tests/bin/make_humann_fixtures.sh:

  stage1 <staging> <out>  reads + database root from HUMAnN 4.0.0a2's bundled
                          demo data (copied out of the pinned image into
                          <staging> by the .sh wrapper)
  stage2 <out>            synthetic KO / EC maps over the gene families the
                          fixture run actually produced (written by the
                          wrapper's real HUMAnN run), plus name-map subsets

Output layout (tests/data/humann/):
  reads/htest_R1.fastq.gz, reads/htest_R2.fastq.gz
                       every 50th read of humann/tests/data/demo.fastq,
                       alternating between the two files (mates are just
                       separate reads to HUMAnN)
  db/chocophlan/       the demo ChocoPhlAn pangenomes (SGB1871/2091/2301)
  db/uniref/           the demo UniRef90 DIAMOND database
  db/utility_mapping/  the demo MetaCyc pathway files + synthetic
                       map_ko_uniref90 / map_level4ec_uniclust90 + name subsets
  db/metaphlan/        EMPTY placeholder mpa_vOct22_CHOCOPhlAnSGB_202403 files:
                       the fixture runs use --bypass-prescreen (no MetaPhlAn
                       database in CI); the placeholders satisfy the module's
                       layout check only
  db/DB_VERSION

The synthetic maps are not biology: they are built so every composer case is
exercised - a UniRef90 family with two KOs, a family with both a KO and an EC,
UniClust90 families with an EC (never a KO, as in the real v4 maps), mapped
and unmapped families, and map members absent from the sample.
"""

import gzip
import shutil
import sys
from pathlib import Path

MPA_INDEX = 'mpa_vOct22_CHOCOPhlAnSGB_202403'


def write_gz(path, text):
    # mtime=0 + no filename: byte-stable across regenerations
    with open(path, 'wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz:
            gz.write(text.encode())


def stage1(staging, out):
    staging, out = Path(staging), Path(out)
    reads = out / 'reads'
    reads.mkdir(parents=True, exist_ok=True)
    with open(staging / 'demo.fastq') as fh:
        lines = fh.read().splitlines()
    records = [lines[i:i + 4] for i in range(0, len(lines), 4)]
    picked = records[::50]
    r1 = [rec for i, rec in enumerate(picked) if i % 2 == 0]
    r2 = [rec for i, rec in enumerate(picked) if i % 2 == 1]
    write_gz(reads / 'htest_R1.fastq.gz', ''.join('\n'.join(r) + '\n' for r in r1))
    write_gz(reads / 'htest_R2.fastq.gz', ''.join('\n'.join(r) + '\n' for r in r2))

    db = out / 'db'
    for sub in ('chocophlan', 'uniref'):
        target = db / sub
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(staging / sub, target)
    util = db / 'utility_mapping'
    util.mkdir(parents=True, exist_ok=True)
    for name in ('metacyc_reactions_level4ec_only.uniref.bz2',
                 'metacyc_pathways_structured_filtered_v24_subreactions'):
        shutil.copyfile(staging / 'utility' / name, util / name)
    mpa = db / 'metaphlan'
    mpa.mkdir(parents=True, exist_ok=True)
    for ext in ('pkl', '1.bt2l', '2.bt2l', '3.bt2l', '4.bt2l', 'rev.1.bt2l', 'rev.2.bt2l'):
        (mpa / f'{MPA_INDEX}.{ext}').write_bytes(b'')
    (db / 'DB_VERSION').write_text(
        'fixture: HUMAnN 4.0.0a2 bundled demo databases + synthetic KO/EC maps '
        '(tests/bin/make_humann_fixtures.py)\n')
    print(f"reads: {len(r1)} + {len(r2)}")


def stage2(out, staging):
    out, staging = Path(out), Path(staging)
    families = []
    with open(out / 'htest_2_genefamilies.tsv') as fh:
        next(fh)
        for line in fh:
            feature = line.split('\t')[0]
            if '|' not in feature and feature != 'READS_UNMAPPED':
                families.append(feature)
    uniref = sorted(f for f in families if f.startswith('UniRef90_'))
    uniclust = sorted(f for f in families if f.startswith('UniClust90_'))
    if len(uniref) < 6 or len(uniclust) < 2:
        sys.exit(f"fixture run produced too few families ({len(uniref)} UniRef90, "
                 f"{len(uniclust)} UniClust90) for the synthetic maps")

    # KO (UniRef90 members only, like the real map_ko_uniref90):
    #   uniref[0] -> K00001 and K00002 (one family, two KOs)
    #   uniref[1], uniref[2] -> K00002 (two families, one KO)
    #   uniref[3] -> K00003; a member absent from the sample rides along
    ko = {
        'K00001': [uniref[0], 'UniRef90_ABSENT1'],
        'K00002': [uniref[0], uniref[1], uniref[2]],
        'K00003': [uniref[3]],
    }
    # EC (UniRef90 and UniClust90 members, like the real map_level4ec_uniclust90):
    #   uniref[0] -> 1.1.1.1 (KO + EC); uniclust[0], uniclust[1] -> 2.7.7.7;
    #   uniref[4] -> 2.7.7.7 and 3.1.1.1 (one family, two ECs)
    ec = {
        '1.1.1.1': [uniref[0]],
        '2.7.7.7': [uniclust[0], uniclust[1], uniref[4], 'UniClust90_ABSENT'],
        '3.1.1.1': [uniref[4]],
    }
    util = out / 'db' / 'utility_mapping'
    write_gz(util / 'map_ko_uniref90.txt.gz',
             ''.join(f"{k}\t" + '\t'.join(v) + '\n' for k, v in ko.items()))
    write_gz(util / 'map_level4ec_uniclust90.txt.gz',
             ''.join(f"{k}\t" + '\t'.join(v) + '\n' for k, v in ec.items()))

    for name, keys in (('map_ko_name.txt.gz', ko), ('map_level4ec_name.txt.gz', ec)):
        with gzip.open(staging / 'utility' / name, 'rt') as fh:
            subset = [line for line in fh if line.split('\t')[0] in keys]
        write_gz(util / name, ''.join(subset))
    print(f"families: {len(uniref)} UniRef90, {len(uniclust)} UniClust90; maps written")


if __name__ == '__main__':
    if sys.argv[1] == 'stage1':
        stage1(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == 'stage2':
        stage2(sys.argv[2], sys.argv[3])
    else:
        sys.exit('usage: make_humann_fixtures.py stage1 <staging> <out> | stage2 <out> <staging>')
