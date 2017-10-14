/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
  CHECDB_Zero Ver.1.0
*/

#include "CHECDB.h"
#include "CHECData.h"

CHECDB_Zero::CHECDB_Zero() : CHECDataBlock(HECDB_ZERO), zero(0) {}

CHECDB_Zero::~CHECDB_Zero() { Clear(); }

void CHECDB_Zero::Clear() { zero = 0; }

void CHECDB_Zero::Write(CHECData *hecd) {
  hecd->WriteHeader("!ZERO");
  hecd->WriteData("F", zero);
}

bool CHECDB_Zero::Read(CHECData *hecd, char *header_line) {
  int rcode[10];
  return hecd->ReadData(rcode, "F", &zero);
}
