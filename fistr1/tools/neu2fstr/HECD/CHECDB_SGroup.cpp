/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
  CHECDB_SGroup Ver.1.0
*/

#include "CHECDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;

CHECDB_SGroup::CHECDB_SGroup() : CHECDataBlock(HECDB_SGROUP), ItemList() {
  name[0] = 0;
}

CHECDB_SGroup::~CHECDB_SGroup() {}

void CHECDB_SGroup::Clear() { ItemList.clear(); }

void CHECDB_SGroup::Write(CHECData *hecd) {
  if (ItemList.size() == 0) return;

  hecd->WriteHeader("!SGROUP", "S", "SGRP", name);
  vector<CItem>::iterator iter;

  for (iter = ItemList.begin(); iter != ItemList.end(); iter++) {
    hecd->WriteData("II", iter->elem, iter->surf);
  }
}

bool CHECDB_SGroup::Read(CHECData *hecd, char *header_line) {
  int rcode[5];

  if (!hecd->ParseHeader(header_line, rcode, "S", "SGRP", name)) return false;

  char line[256];
  const int max_id_n = 100;
  int id[max_id_n];
  int i, n, n2;

  while (1) {
    if (!hecd->ReadLine(line)) break;

    if (line[0] == '!') {
      hecd->PushReadLine(line);
      break;
    }

    n = hecd->ParseIntDataArray(line, id);

    if (n < 0) return false;

    n2 = n / 2;

    if (n2 % 2 != 0) return false;

    for (i = 0; i < n2; i++) {
      ItemList.push_back(CItem(id[2 * i], id[2 * i + 1]));
    }
  }

  return true;
}
