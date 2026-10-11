#!/usr/bin/env python3
"""Delegating wrapper around eggNOG-mapper's emapper.py (design doc Q23).

eggNOG-mapper 3.0.0-beta6 parses the OG memberships of a seed ortholog with
AnnotationEngine._parse_ogs_string, which splits each entry on every '|'.
eggNOG 7 has OGs whose family name itself contains '|' (2.3 % of OG names,
e.g. 'ABC_tran|TL31Y9@131567|A-1'); the method truncates those keys
('ABC_tran|TL31Y9@131567', level '-'), so their metadata lookup fails and
the gene loses or mis-assigns its COG category (upstream issue #620, open).

This wrapper replaces the method with the candidate fix posted on that
issue — split the family from '@taxid' first, then read the pipe-delimited
fields after the taxid; anything else keeps the original split — and then
executes the stock CLI unchanged, forwarding all arguments. Annotation
workers are forked from this process, so they inherit the patched class.

The patch applies only to the exact beta6 method: if the installed source
differs (an upstream fix or any other change), the wrapper stops with an
error instead of patching blindly. Delete it at the eggNOG-mapper re-pin
(design doc Q11) once upstream has fixed #620.
"""

import hashlib
import inspect
import runpy
import shutil
import sys

# sha256 of inspect.getsource(AnnotationEngine._parse_ogs_string) as shipped
# in eggNOG-mapper v3.0.0-beta6 (tag commit b3757a6)
BETA6_METHOD_SHA256 = '7b2c4b0a8419845d6807be74b7a5efae79966c02ce4e9d5da064d002b93255d1'


def parse_ogs_string(self, ogs_string):
    """Fixed copy of AnnotationEngine._parse_ogs_string (upstream #620).

    Format: "family@taxid|clade|taxon_name,family@taxid|clade|taxon_name,..."
    where the family name itself may contain '|'. Returns a list of
    (og_name, level) with og_name = "family@taxid|clade" and level = taxid.
    """
    result = []
    for og in ogs_string.split(","):
        og = og.strip()
        if not og:
            continue
        # Family names can contain "|"; parse the fields after "@taxid"
        family, at, suffix = og.partition("@")
        fields = suffix.split("|")
        if at and len(fields) >= 2 and fields[0].isdigit() and fields[1]:
            parts = [family + at + fields[0], fields[1]]
        else:
            # Original behaviour for entries that do not fit the format
            parts = og.split("|")
        if len(parts) >= 2:
            og_name = "|".join(parts[:2])
            if "@" in parts[0]:
                level = parts[0].split("@")[1]
            else:
                level = "-"
            result.append((og_name, level))
    return result


def patch(engine_cls):
    """Replace engine_cls._parse_ogs_string after checking it is the beta6 one."""
    try:
        installed = inspect.getsource(engine_cls._parse_ogs_string)
    except (OSError, TypeError):
        installed = ''   # no source to verify: treat as not the beta6 method
    digest = hashlib.sha256(installed.encode()).hexdigest()
    if digest != BETA6_METHOD_SHA256:
        sys.exit(
            'emapper_ogs_fix.py: the installed eggNOG-mapper '
            '_parse_ogs_string is not the v3.0.0-beta6 method this patch '
            f'targets (sha256 {digest}). Check whether upstream issue #620 '
            'is fixed in this version, then update or delete '
            'bin/emapper_ogs_fix.py (design doc Q23).')
    engine_cls._parse_ogs_string = parse_ogs_string


def main():
    from eggnogmapper.annotator.e7.annotate import AnnotationEngine

    patch(AnnotationEngine)

    cli = shutil.which('emapper.py')
    if cli is None:
        sys.exit('emapper.py not found on PATH — is this running inside the '
                 'eggNOG-mapper container?')
    sys.argv[0] = cli
    runpy.run_path(cli, run_name='__main__')


if __name__ == '__main__':
    main()
