/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
  CFSTRDB_Reftemp Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;

CFSTRDB_Reftemp::CFSTRDB_Reftemp() : CFSTRDataBlock(FSTRDB_REFTEMP), value(0) {}

CFSTRDB_Reftemp::~CFSTRDB_Reftemp() { Clear(); }

void CFSTRDB_Reftemp::Clear() { value = 0; }

void CFSTRDB_Reftemp::Write(CHECData *hecd) {
  hecd->WriteHeader("!REFTEMP");
  hecd->WriteData("F", value);
}

bool CFSTRDB_Reftemp::Read(CHECData *hecd, char *header_line) {
  int rcode[5];
  return hecd->ReadData(rcode, "F", &value);
}
