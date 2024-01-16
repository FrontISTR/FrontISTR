#!/bin/bash

#> git clone
git submodule update --init --recursive
cd hecmw1/monolis

#> install monoils_utils
cd submodule/monolis_utils
make lib

#> install gedatsu
cd ../gedatsu
make FLAGS=SUBMODULE lib

#> install ggtools
cd ../ggtools
make FLAGS=SUBMODULE lib

#> install monoils
cd ../..
make lib
