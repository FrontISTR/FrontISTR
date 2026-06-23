/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
  CFSTRDB_SRadiate Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;

CFSTRDB_SRadiate::CFSTRDB_SRadiate()
    : CFSTRDataBlock(FSTRDB_SRADIATE), ItemList() {
  amp1[0] = 0;
  amp2[0] = 0;
}

CFSTRDB_SRadiate::~CFSTRDB_SRadiate() { Clear(); }

void CFSTRDB_SRadiate::Clear() {
  ItemList.clear();
  amp1[0] = 0;
  amp2[0] = 0;
}

void CFSTRDB_SRadiate::Write(CHECData *hecd) {
  char buff[256];

  if (ItemList.size() == 0) return;

  size_t len = snprintf(buff, sizeof(buff), "!SRADIATE");

  if (amp1[0] != 0) {
    len += snprintf(buff + len, sizeof(buff) - len, ",AMP1=%s", amp1);
  }

  if (amp2[0] != 0) {
    snprintf(buff + len, sizeof(buff) - len, ",AMP2=%s", amp2);
  }

  hecd->WriteHeader(buff);
  vector<CItem>::iterator iter;

  for (iter = ItemList.begin(); iter != ItemList.end(); iter++) {
    hecd->WriteData("SFF", iter->sgrp, iter->value, iter->sink);
  }
}

bool CFSTRDB_SRadiate::Read(CHECData *hecd, char *header_line) {
  int rcode[10];
  amp1[0] = 0;
  amp2[0] = 0;

  if (!hecd->ParseHeader(header_line, rcode, "SS", "AMP1", amp1, "AMP2", amp2))
    return false;

  while (1) {
    CItem item;
    bool fg = hecd->ReadData(rcode, "SFF", item.sgrp, &item.value, &item.sink);

    if (!fg) break;

    ItemList.push_back(item);
  }

  return true;
}
