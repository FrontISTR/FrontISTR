/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
  CHECDB_Amplitude Ver.1.0
*/

#include "CHECDB.h"
#include "CHECData.h"

using namespace std;

CHECDB_Amplitude::CHECDB_Amplitude()
    : CHECDataBlock(HECDB_AMPLITUDE), ItemList() {
  name[0] = 0;
  strcpy(definition, "TABULAR");
  strcpy(time, "STEP TIME");
  strcpy(value, "RLATIVE");
}

CHECDB_Amplitude::~CHECDB_Amplitude() {}

void CHECDB_Amplitude::Clear() { ItemList.clear(); }

void CHECDB_Amplitude::Write(CHECData *hecd) {
  if (ItemList.size() == 0) return;

  hecd->WriteHeader("!AMPLITUDE", "SSSS", "NAME", name, "DEFINITION",
                    definition, "TIME", time, "VALUE", value);
  vector<CItem>::iterator iter;

  for (iter = ItemList.begin(); iter != ItemList.end(); iter++) {
    hecd->WriteData("FF", iter->val, iter->t);
  }
}

bool CHECDB_Amplitude::Read(CHECData *hecd, char *header_line) {
  int rcode[10];

  if (!hecd->ParseHeader(header_line, rcode, "SSSS", "NAME", name, "DEFINITION",
                         definition, "TIME", time, "VALUE", value))
    return false;

  double val, t;

  while (hecd->ReadData(rcode, "FF", &val, &t)) {
    ItemList.push_back(CItem(val, t));
  }

  return true;
}
