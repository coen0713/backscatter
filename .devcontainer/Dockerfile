FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential g++-14 clang-18 clang-tidy-18 lldb-18 \
      cmake ninja-build git ca-certificates \
      libgdal-dev gdal-bin \
      python3 python3-pip python3-venv \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir numpy matplotlib requests clang-format==19.1.7 nanobind
ENV PATH="/opt/venv/bin:${PATH}"

WORKDIR /workspace
