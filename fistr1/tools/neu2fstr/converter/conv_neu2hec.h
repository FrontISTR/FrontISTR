/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
        conv_neu2hec ver.1.0
*/

#ifndef conv_neu2hecH
#define conv_neu2hecH

#include "CNFData.h"
#include "CNFDataBlock.h"
#include "CHECData.h"
#include "CHECDB.h"

// solution
enum { sol_static = 0, sol_heat, sol_eigen };

void conv_neu2hec(CNFData& neu, CHECData& hec, int solution = 0);

#endif
