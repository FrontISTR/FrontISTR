/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "hecmw_config.h"
#include "hecmw_msg.h"
#include "hecmw_util.h"

static struct hecmw_msgent msg_unknown = {-1, "HEC-MW-UNKNOWN",
                                          "Unkown message No"};

static struct hecmw_msgent msg_err = {-1, "HEC-MW-SYSERR", NULL};

static char msg_buf[HECMW_MSG_LEN + 1];

static struct hecmw_msgent *get_msgent(int msgno) {
  int i;
  struct hecmw_msgent *p;

  i = 0;
  for (p = &hecmw_msg_table[i]; p->msgno != -1; p = &hecmw_msg_table[i++]) {
    if (msgno == p->msgno) return p;
  }

  return NULL; /* not found */
}

char *HECMW_strmsg(int msgno) {
  struct hecmw_msgent *p;

  p = NULL;
  if (msgno < HECMW_MSGNO_BASE) {
    /* system error */
    p      = &msg_err;
    p->msg = strerror(msgno);
  }
  if (p == NULL) p = get_msgent(msgno);
  if (p == NULL) p = &msg_unknown;
  sprintf(msg_buf, "%s: %s(%d)", p->msgno_str, p->msg, msgno);
  return msg_buf;
}

int HECMW_is_syserr(int msgno) { return msgno < HECMW_MSGNO_BASE; }

void hecmw_strmsg_if(int *msgno, char *dst, int dstlen) {
  const char *p;

  if (dst == NULL) return;
  if (dstlen < 0) return;

  p = HECMW_strmsg(*msgno);
  HECMW_strcpy_c2f(p, dst, dstlen);
}

void hecmw_strmsg_if_(int *msgno, char *dst, int dstlen) {
  hecmw_strmsg_if(msgno, dst, dstlen);
}

void hecmw_strmsg_if__(int *msgno, char *dst, int dstlen) {
  hecmw_strmsg_if(msgno, dst, dstlen);
}

void HECMW_STRMSG_IF(int *msgno, char *dst, int dstlen) {
  hecmw_strmsg_if(msgno, dst, dstlen);
}
