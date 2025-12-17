#!/bin/bash

source ~/.bashrc
set -x
set -e
umask 022
# exit 1 if error

## Spack install spec for desired compiler
#SpackCompiler="gcc@12.2.0"

spack install gcc@12.2.0

#spack compiler find --scope system
#spack external find --scope system
#spack install -j32 --fail-fast $SpackCompiler --show-log-on-error target=x86_64 platform=linux os=amzn2
#spack load $SpackCompiler
#spack compiler find --scope system
