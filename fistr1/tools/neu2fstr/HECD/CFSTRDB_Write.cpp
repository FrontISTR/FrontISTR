/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
  CFSTRDB_Write Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;

CFSTRDB_Write::CFSTRDB_Write()
    : CFSTRDataBlock(FSTRDB_WRITE), result(0), visual(0) {}

CFSTRDB_Write::~CFSTRDB_Write() { Clear(); }

void CFSTRDB_Write::Clear() { result = visual = 0; }

void CFSTRDB_Write::Write(CHECData *hecd) {
  char header_s[256];
  strcpy(header_s, "!WRITE");

  if (result) strcat(header_s, ",RESULT");

  if (visual) strcat(header_s, ",VISUAL");

  hecd->WriteHeader(header_s);
}

bool CFSTRDB_Write::Read(CHECData *hecd, char *header_line) {
  int rcode[10];

  if (!hecd->ParseHeader(header_line, rcode, "EE", "RESULT", &result, "VISUAL",
                         &visual))
    return false;

  return true;
}
