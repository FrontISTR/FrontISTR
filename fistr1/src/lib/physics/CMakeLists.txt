###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr
	PRIVATE
    material.f90
    mechgauss.f90
    ElasticLinear.f90
    Hyperelastic.f90
    Elastoplastic.f90
    Viscoelastic.f90
    creep.f90
    calMatMatrix.f90
)
