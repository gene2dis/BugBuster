# eggNOG-mapper 3.0.0-beta6 for docker/podman (design doc Q11 interim, 2026-10-10).
#
# Upstream publishes v3 only as an Apptainer image
# (https://data.cgmlab.org/eggnog-mapper/emapper-3.0/eggnog-mapper-3.0.0-beta6.sif),
# which the pipeline keeps using under singularity/apptainer. This is a port
# of upstream's own recipe (apptainer/eggnog-mapper.def at tag v3.0.0-beta6)
# pinned to what that .sif contains, so both engines run the same code:
#   - base python:3.12.13-slim-bookworm (the .sif: Debian 12, Python 3.12.13)
#   - DIAMOND 2.1.24, MMseqs2 commit 17b688d and Prodigal 2.6.3: the exact
#     binaries in the .sif (sha256-verified below; upstream's recipe takes
#     MMseqs2 'latest', so the archived build of the .sif's commit is used)
#   - HMMER 3.4 compiled from source, as upstream does
#   - Python packages at the versions installed in the .sif
#   - eggNOG-mapper from the v3.0.0-beta6 tag commit
# One deliberate addition: procps (`ps`), which Nextflow needs for task
# metrics and the .sif lacks (there, bin/ps stands in).
#
# The upstream #620 fix is NOT baked in: the pipeline applies it at run time
# under both engines (bin/emapper_ogs_fix.py, design doc Q23).
#
# Built and pushed to ghcr.io/gene2dis/bugbuster-eggnog-mapper by
# .github/workflows/eggnog-image.yml. Retire this image when eggNOG-mapper
# v3.0.0 final reaches bioconda (design doc Q11).

FROM python:3.12.13-slim-bookworm@sha256:4766d8b510c428e595d74b9cc5bbb2fae8e26316fffb4adc89908d79aacd58a2

ARG EMAPPER_COMMIT=b3757a6d226047527a729546e58ff530d76a5d7d
ARG DIAMOND_VERSION=2.1.24
ARG DIAMOND_SHA256=5d04a9ecbbceea5d1a6496070155012c2538788def552f85e0d65095a532a5f5
ARG MMSEQS_COMMIT=17b688d21dda57fc5f5b7286ecba7ec003d4717f
ARG MMSEQS_SHA256=f35d2a3642987da84457a644be25090b764ec8ffe1f3f0f79c5ec64bc21709d2
ARG HMMER_VERSION=3.4
ARG PRODIGAL_VERSION=2.6.3
ARG PRODIGAL_SHA256=80a3d7acb08ba9f63dde6b6ef4e95ec356128ac355f36960ea722f381150f246

LABEL org.opencontainers.image.title="eggNOG-mapper" \
      org.opencontainers.image.version="3.0.0-beta6" \
      org.opencontainers.image.description="eggNOG-mapper 3.0.0-beta6 (upstream recipe, .sif-pinned tool versions, plus procps) for the BugBuster pipeline" \
      org.opencontainers.image.source="https://github.com/gene2dis/BugBuster" \
      org.opencontainers.image.licenses="AGPL-3.0"

ENV PATH=/opt/diamond:/opt/mmseqs/bin:/opt/hmmer/bin:/opt/prodigal:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Base system, as upstream, plus procps
RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
        ca-certificates curl gcc g++ make zlib1g-dev libgomp1 procps \
    && rm -rf /var/lib/apt/lists/*

# DIAMOND, MMseqs2 and Prodigal: the binaries shipped in the .sif
RUN mkdir -p /opt/diamond /opt/prodigal \
    && curl -fsSL -o /tmp/diamond.tar.gz \
        "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz" \
    && tar -xzf /tmp/diamond.tar.gz -C /opt/diamond \
    && echo "${DIAMOND_SHA256}  /opt/diamond/diamond" | sha256sum -c - \
    && curl -fsSL -o /tmp/mmseqs.tar.gz \
        "https://mmseqs.com/archive/${MMSEQS_COMMIT}/mmseqs-linux-avx2.tar.gz" \
    && tar -xzf /tmp/mmseqs.tar.gz -C /opt \
    && echo "${MMSEQS_SHA256}  /opt/mmseqs/bin/mmseqs" | sha256sum -c - \
    && curl -fsSL -o /opt/prodigal/prodigal \
        "https://github.com/hyattpd/Prodigal/releases/download/v${PRODIGAL_VERSION}/prodigal.linux" \
    && echo "${PRODIGAL_SHA256}  /opt/prodigal/prodigal" | sha256sum -c - \
    && chmod +x /opt/diamond/diamond /opt/prodigal/prodigal \
    && rm -f /tmp/diamond.tar.gz /tmp/mmseqs.tar.gz

# HMMER, compiled from source as upstream does
RUN mkdir -p /tmp/hmmer-build \
    && curl -fsSL "http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz" \
        | tar -xz -C /tmp/hmmer-build --strip-components=1 \
    && cd /tmp/hmmer-build \
    && ./configure --prefix=/opt/hmmer > /dev/null \
    && make -j"$(nproc)" > /dev/null \
    && make install > /dev/null \
    && cd / && rm -rf /tmp/hmmer-build

# Python packages at the .sif's versions, then eggNOG-mapper from the tag commit
RUN python3.12 -m pip install --no-cache-dir \
        cython==3.2.9 biopython==1.88 numpy==2.5.2 psutil==7.2.2 xlsxwriter==3.2.9 \
    && curl -fsSL -o /tmp/emapper-src.tar.gz \
        "https://github.com/eggnogdb/eggnog-mapper/archive/${EMAPPER_COMMIT}.tar.gz" \
    && mkdir -p /opt/emapper-src \
    && tar -xzf /tmp/emapper-src.tar.gz -C /opt/emapper-src --strip-components=1 \
    && python3.12 -m pip install --no-cache-dir /opt/emapper-src \
    && rm -rf /opt/emapper-src /tmp/emapper-src.tar.gz

# Upstream's %test smoke checks
RUN emapper.py --version \
    && diamond --version \
    && mmseqs version \
    && hmmsearch -h | head -2 \
    && prodigal -v 2>&1 | head -2 \
    && python3.12 -c "import eggnogmapper.annotator._codec, eggnogmapper.annotator._collect_inner; print('Cython extensions loaded OK')" \
    && ps -e -o pid= -o ppid= > /dev/null

CMD ["emapper.py", "--help"]
