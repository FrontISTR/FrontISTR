/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
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
  size_t len = snprintf(header_s, sizeof(header_s), "!WRITE");

  if (result)
    len += snprintf(header_s + len, sizeof(header_s) - len, ",RESULT");

  if (visual) snprintf(header_s + len, sizeof(header_s) - len, ",VISUAL");

  hecd->WriteHeader(header_s);
}

bool CFSTRDB_Write::Read(CHECData *hecd, char *header_line) {
  int rcode[10];

  if (!hecd->ParseHeader(header_line, rcode, "EE", "RESULT", &result, "VISUAL",
                         &visual))
    return false;

  return true;
}
