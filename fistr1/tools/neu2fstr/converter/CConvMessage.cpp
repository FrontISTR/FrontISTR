/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/* CConvMessage class Ver.1.0 */

#include "CConvMessage.h"

const char ERROR_MSG[][80] = {"No error",
                              "Unknown error",
                              "Coordinate error",
                              "Not supported element",
                              "Invalid element property",
                              "Not supported property of element"};

char CConvMessage::msg[256] = "";

CConvMessage::CConvMessage(int No, const char *op_msg, ...) : no(No) {
  if (op_msg[0] == 0) {
    option_msg[0] = 0;
    return;
  }

  va_list va;
  va_start(va, op_msg);
  vsprintf(option_msg, op_msg, va);
  va_end(va);
}

const char *CConvMessage::Msg() {
  if (option_msg[0] != 0) {
    snprintf(msg, sizeof(msg), "##Error: %s : %s", ERROR_MSG[no], option_msg);

  } else {
    snprintf(msg, sizeof(msg), "##Error: %s", ERROR_MSG[no]);
  }

  return msg;
}
