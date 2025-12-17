#!/bin/bash

source ~/.bashrc
set -x
set -e
umask 022 # exit 1 if error

## Spack install spec for desired compiler
SpackCompiler="gcc@12.2.0"

spack compiler find
spack install -j32 --fail-fast $SpackCompiler --show-log-on-error target=x86_64 platform=linux os=ubuntu24.04
spack load $SpackCompiler

# Add compiler to ./spack/linux/compilers.yaml
spack compiler find
