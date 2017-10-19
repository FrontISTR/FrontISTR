/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
  CFSTRDB_Eigen Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;

CFSTRDB_Eigen::CFSTRDB_Eigen()
    : CFSTRDataBlock(FSTRDB_EIGEN), nset(5), lcztol(1e-8), lczmax(60) {}

CFSTRDB_Eigen::~CFSTRDB_Eigen() { Clear(); }

void CFSTRDB_Eigen::Clear() {
  nset   = 5;
  lcztol = 1e-8;
  lczmax = 60;
}

void CFSTRDB_Eigen::Write(CHECData *hecd) {
  hecd->WriteHeader("!EIGEN");
  hecd->WriteData("IFI", nset, lcztol, lczmax);
}

bool CFSTRDB_Eigen::Read(CHECData *hecd, char *header_line) {
  int rcode[5];
  return hecd->ReadData(rcode, "IFI", &nset, &lcztol, &lczmax);
}
