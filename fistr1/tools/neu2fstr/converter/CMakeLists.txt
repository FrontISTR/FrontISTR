###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(neu2fstr
  PRIVATE
    cconv_mat.cpp
    CConvMessage.cpp
    conv_neu2hec.cpp
    conv_neu2fstr_static.cpp
    conv_neu2fstr_dynamic.cpp
    conv_neu2fstr_heat.cpp
  PUBLIC
    CConvMessage.h
    cconv_mat.h
    conv_neu2fstr_dynamic.h
    conv_neu2fstr_heat.h
    conv_neu2fstr_static.h
    conv_neu2hec.h
    conv_util.h
)

