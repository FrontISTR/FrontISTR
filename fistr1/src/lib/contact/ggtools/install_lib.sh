#!/bin/bash

git submodule update --init --recursive

BASE_DIR=$(pwd)

#> monolis_utils
cd submodule/monolis_utils
make
#make FLAGS=INTEL
