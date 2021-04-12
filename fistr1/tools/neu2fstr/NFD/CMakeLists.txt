###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
target_sources(neu2fstr
  PRIVATE
    CNFData.cpp
    CNFDB_100.cpp
    CNFDB_402.cpp
    CNFDB_403.cpp
    CNFDB_404.cpp
    CNFDB_405.cpp
    CNFDB_408.cpp
    CNFDB_506.cpp
    CNFDB_507.cpp
    CNFDB_601.cpp
    CNFMessage.cpp
  PUBLIC
    CNFDB_100.h
    CNFDB_402.h
    CNFDB_403.h
    CNFDB_404.h
    CNFDB_405.h
    CNFDB_408.h
    CNFDB_506.h
    CNFDB_507.h
    CNFDB_601.h
    CNFData.h
    CNFDataBlock.h
    CNFMessage.h
)

