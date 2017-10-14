/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/* CConvMessage class Ver. 3.6 */

#ifndef CConvMessageH
#define CConvMessageH

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

enum {
  CONV_NO_ERROR = 0,
  CONV_UNKNOWN_ERROR,
  CONV_COORDINATE_ERROR,
  CONV_NO_SUPPORTED_ELEMENT,
  CONV_INVALID_ELEMENT_PROPERTY,
  CONV_NO_SUPPORTED_ELEM_PROP
};

class CConvMessage {
 protected:
  static char msg[256];

 public:
  char option_msg[256];
  int no;
  CConvMessage(int No = 0, const char* op_msg = "", ...);
  virtual const char* Msg();
};

#endif
