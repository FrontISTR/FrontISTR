/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_MSG_INCLUDED
#define HECMW_MSG_INCLUDED

#include "hecmw_msgno.h"

struct hecmw_msgent {
  int msgno;
  char *msgno_str;
  char *msg;
};

extern struct hecmw_msgent hecmw_msg_table[];

extern char *HECMW_strmsg(int msgno);

extern int HECMW_is_syserr(int msgno);

#endif
