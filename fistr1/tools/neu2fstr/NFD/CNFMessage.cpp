/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
  CNFMessage ver.1.0
*/

#include <stdio.h>
#include "CNFMessage.h"

static char ERROR_MSG[][80] = {
    "No error",        "Unknown error",      "Cannot open NEU file",
    "Read data error", "Data line required", "Data block required",
    "Invalid token",   "Item required",      "A record is required"};

static char WARNING_MSG[][80] = {"Non supported data block"};

char CNFMessage::msg[256] = "";

const char *CNFError::Msg() {
  size_t len = snprintf(msg, sizeof(msg), "##Error");

  if (line >= 0) {
    len += snprintf(msg + len, sizeof(msg) - len, "(line:%d", line);

    if (column > 0) {
      len += snprintf(msg + len, sizeof(msg) - len, ",col:%d", column);
    }

    len += snprintf(msg + len, sizeof(msg) - len, ")");
  }

  snprintf(msg + len, sizeof(msg) - len, ": %s%s", ERROR_MSG[no], option_msg);
  return msg;
}

const char *CNFWarning::Msg() {
  size_t len = snprintf(msg, sizeof(msg), "##Warning");

  if (line >= 0) {
    len += snprintf(msg + len, sizeof(msg) - len, "(line:%d", line);

    if (column > 0) {
      len += snprintf(msg + len, sizeof(msg) - len, ",col:%d", column);
    }

    len += snprintf(msg + len, sizeof(msg) - len, ")");
  }

  snprintf(msg + len, sizeof(msg) - len, ": %s%s", WARNING_MSG[no], option_msg);
  return msg;
}
