# Docker equivalent of Nanopore_baseContainer.sing, for running ViralFlow's
# NANOPORE mode where Singularity/Apptainer is unavailable - CI runners and
# developer machines, notably Apple Silicon.
#
# Build and use with the `docker` profile:
#
#   cd vfnext/containers
#   docker build -f nanopore_base.Dockerfile -t viralflow/nanopore-base:2.0.0a1 .
#   nextflow run ../main.nf -profile docker --mode NANOPORE ...
#
# The version pins below MUST match Nanopore_baseContainer.sing. They are two
# recipes for one environment, so tests/test_container_recipes.py fails the
# build if they drift apart. See that file before changing a pin here.
FROM ubuntu:24.04

ARG HTSLIB_VERSION=1.21
ARG MINIMAP2_TAG=v2.28
ARG NETWORKX_VERSION=3.6.1
# Pinned to commits rather than tags, on purpose - see Nanopore_baseContainer.sing
# for the reasoning (Porechop_ABI master is ahead of its last tag; bamUtil's last
# tag clones libStatGen over the retired git:// protocol).
ARG PORECHOP_ABI_COMMIT=0bc9f17f31d4ec1dcbab4796871cc09324cc143b
ARG BAMUTIL_COMMIT=017721cc07948558395e4934ec10d0f91407c5eb

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/app/minimap2/:$PATH

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        curl \
        wget \
        less \
        ca-certificates \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-dev \
        gcc \
        # automake provides aclocal, which autoreconf needs. The .sing gets it
        # implicitly as an apt recommendation of autoconf; naming it here keeps
        # the build independent of how Ubuntu tunes its recommends.
        make autoconf automake \
        libbz2-dev liblzma-dev libncurses5-dev \
        libcurl4-openssl-dev libssl-dev zlib1g-dev \
        libgsl-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# minimap2. Upstream's plain `make` targets x86 SSE; on aarch64 it needs the
# NEON flags, which is why this differs from the .sing (that recipe is only ever
# built on x86_64).
RUN git clone --depth 1 --branch ${MINIMAP2_TAG} https://github.com/lh3/minimap2 \
    && cd minimap2 \
    && if [ "$(uname -m)" = "aarch64" ]; then make arm_neon=1 aarch64=1; else make; fi

# htslib
RUN wget -O htslib.tar.bz2 "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
    && tar -xf htslib.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && autoheader && autoreconf -i \
    && ./configure --prefix=/usr/local/ \
    && make && make install \
    # Refresh the linker cache. bcftools is built from the GitHub source
    # archive, which does not bundle htslib, so it links dynamically against
    # /usr/local/lib/libhts.so. Without this it builds fine and then fails at
    # runtime with "libhts.so.3: cannot open shared object file".
    && ldconfig \
    && cd /app && rm -rf htslib.tar.bz2 htslib-${HTSLIB_VERSION}

# bcftools
RUN wget -O bcftools.tar.gz "https://github.com/samtools/bcftools/archive/${HTSLIB_VERSION}.tar.gz" \
    && tar xf bcftools.tar.gz \
    && cd bcftools-${HTSLIB_VERSION} \
    && autoheader && autoconf \
    && ./configure --enable-libgsl \
    && make \
    && mv bcftools /usr/local/bin \
    && cd /app && rm -rf bcftools.tar.gz bcftools-${HTSLIB_VERSION}

# samtools
RUN wget -O samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/${HTSLIB_VERSION}/samtools-${HTSLIB_VERSION}.tar.bz2" \
    && tar xf samtools.tar.bz2 \
    && cd samtools-${HTSLIB_VERSION} \
    && autoheader && autoconf -Wno-syntax \
    && ./configure --prefix /usr/local \
    && make && make install \
    && cd /app && rm -rf samtools.tar.bz2 samtools-${HTSLIB_VERSION}

# Porechop_ABI
RUN pip install --break-system-packages "networkx==${NETWORKX_VERSION}" \
    && git clone https://github.com/bonsai-team/Porechop_ABI.git \
    && cd Porechop_ABI \
    && git checkout ${PORECHOP_ABI_COMMIT} \
    && python3 setup.py install

# bamUtil. `make cloneLib` fetches libStatGen from its own default branch, so
# that dependency remains unpinned here exactly as in the .sing.
RUN git clone https://github.com/statgen/bamUtil.git \
    && cd bamUtil \
    && git checkout ${BAMUTIL_COMMIT} \
    && make cloneLib \
    && make \
    && make install

# Mirrors the %test section of Nanopore_baseContainer.sing: fail the build here
# rather than in the middle of a pipeline run.
RUN python3 --version \
    && minimap2 --version \
    && samtools --version | head -n 1 \
    && bcftools --version | head -n 1 \
    && porechop_abi --version \
    && bam help > /dev/null 2>&1 || true

WORKDIR /
