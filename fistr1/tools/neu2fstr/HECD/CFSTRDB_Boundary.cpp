/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
  CFSTRDB_Boundary Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;

CFSTRDB_Boundary::CFSTRDB_Boundary()
    : CFSTRDataBlock(FSTRDB_BOUNDARY), ItemList() {}

CFSTRDB_Boundary::~CFSTRDB_Boundary() { Clear(); }

void CFSTRDB_Boundary::Clear() { ItemList.clear(); }

void CFSTRDB_Boundary::Write(CHECData *hecd) {
  if (ItemList.size() == 0) return;

  hecd->WriteHeader("!BOUNDARY");
  vector<CItem>::iterator iter;

  for (iter = ItemList.begin(); iter != ItemList.end(); iter++) {
    hecd->WriteData("SIIF", iter->ngrp, iter->dof_ids, iter->dof_ide,
                    iter->value);
  }
}

bool CFSTRDB_Boundary::Read(CHECData *hecd, char *header_line) {
  int rcode[10];

  while (1) {
    CItem item;

    if (!hecd->ReadData(rcode, "SIIF", &item.ngrp, &item.dof_ids, &item.dof_ide,
                        &item.value))
      break;

    ItemList.push_back(item);
  }

  return true;
}
