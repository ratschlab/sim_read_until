# Write a Dockerfile that starts from the latest ubuntu, installs Python3.10 and then installs this package
FROM ubuntu:latest as builder

ARG package_version=unknown_version_docker
SHELL ["/bin/bash", "-c"]

# git for hatch build version, libbz2-dev required for NanoSim
RUN apt-get update && apt-get install -y python3.10 python3-pip git python3.10-venv curl libbz2-dev \
    && mkdir ~/ont_project_all

# # copy single file to avoid rebuilding everything during testing
# # ~/ont_project_all/ont_project/usecases/install_usecase_deps.sh
# COPY usecases/install_usecase_deps.sh install_usecase_deps.sh
# todo2: git clone --depth 1 git@github.com:ratschlab/sim_read_until.git ont_project
# move down while testing to avoid rebuilding everything
# docker does not expand ~, so we need to use /root

# copy installation script to avoid rebuilding mamba environment during testing
RUN mkdir -p /root/install_deps/external/ont_nanosim/ && \
    cd ~/ont_project_all && python3.10 -m venv ont_project_venv \
    && source ~/ont_project_all/ont_project_venv/bin/activate \
    && pip install --upgrade pip
COPY usecases/install_usecase_deps.sh /root/install_deps
COPY external/ont_nanosim/requirements.txt /root/install_deps/external/ont_nanosim/requirements.txt

# install mamba and NanoSim/readfish (conda is unacceptably slow)
# from https://docs.conda.io/en/latest/miniconda.html#linux-installers
# curl -OL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#     && bash Miniconda3-latest-Linux-x86_64.sh -b \
#     && echo 'export PATH="/root/miniconda3/bin:$PATH"' >> ~/.bashrc \
# see https://github.com/conda-forge/miniforge#mambaforge
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh" \
    && bash Mambaforge-$(uname)-$(uname -m).sh -b -p ~/conda \
    && source "${HOME}/conda/etc/profile.d/conda.sh" \
    && source "${HOME}/conda/etc/profile.d/mamba.sh" \
    # skip ReadFish installation
    && (cd /root/install_deps && bash install_usecase_deps.sh) \
    && echo 'export PATH="/root/ont_project_all/tools/bin:$PATH"' >> ~/.bashrc \
    && rm -rf /root/install_deps

COPY . /root/ont_project_all/ont_project

FROM builder as prod
RUN apt-get clean
# no .git directory is copied, so we replace by version based on timestamp
# don't include dots, to make it a valid package version, as described here: https://peps.python.org/pep-0440/, don't use underscores
RUN cd ~/ont_project_all \
    && sed -i'.backup' "s@unknown_version@0.0.1.dev+dockerbuild$(date +%YY%mM%dD%HH%MM%SS)@" ont_project/pyproject.toml \
    && source ~/ont_project_all/ont_project_venv/bin/activate \
    && pip install ~/ont_project_all/ont_project[readfish] \
    && mv ont_project/pyproject.toml.backup ont_project/pyproject.toml \
    && python -c "import simreaduntil" \
    && pip freeze

FROM builder as dev
RUN apt-get install -y liblzma-dev libcurl4-openssl-dev

# no .git directory is copied, so we replace by version based on timestamp
# don't include dots, to make it a valid package version, as described here: https://peps.python.org/pep-0440/, don't use underscores
RUN cd ~/ont_project_all \
    && sed -i'.backup' "s@unknown_version@0.0.1.dev+dockerbuild$(date +%YY%mM%dD%HH%MM%SS)@" ont_project/pyproject.toml \
    && source ~/ont_project_all/ont_project_venv/bin/activate \
    # not installing test target as it requires pytest-sugar which does not work in github runners
    && pip install -e ~/ont_project_all/ont_project[readfish,dev] \
    && mv ont_project/pyproject.toml.backup ont_project/pyproject.toml \
    && python -c "import simreaduntil" \
    && pip freeze

# cd ./ont_project
# docker build --target prod -t ont_simulator:v1 .
# docker build --target dev -t ont_simulator:v1_dev .
# docker build --target dev -t ont_simulator:v1_dev . --progress=plain
# docker run --rm -it ont_simulator:v1
# test Docker image
# docker run --rm -it ont_simulator:v1 -- tox run -- -k test_utils