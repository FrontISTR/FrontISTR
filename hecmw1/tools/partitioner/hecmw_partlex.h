/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_PARTLEX_INCLUDED
#define HECMW_PARTLEX_INCLUDED

#include <stdio.h>

#define HECMW_PARTLEX_NL 1001

#define HECMW_PARTLEX_INT 1002

#define HECMW_PARTLEX_DOUBLE 1003

#define HECMW_PARTLEX_NAME 1004

#define HECMW_PARTLEX_FILENAME 1005

#define HECMW_PARTLEX_H_PARTITION 2001

#define HECMW_PARTLEX_K_TYPE 3001

#define HECMW_PARTLEX_K_METHOD 3011

#define HECMW_PARTLEX_K_DOMAIN 3021

#define HECMW_PARTLEX_K_DEPTH 3031

#define HECMW_PARTLEX_K_UCD 3041

#define HECMW_PARTLEX_K_CONTACT 3051

#define HECMW_PARTLEX_K_PART 3061

#define HECMW_PARTLEX_V_NODE_BASED 3002

#define HECMW_PARTLEX_V_ELEMENT_BASED 3003

#define HECMW_PARTLEX_V_RCB 3012

#define HECMW_PARTLEX_V_KMETIS 3013

#define HECMW_PARTLEX_V_PMETIS 3014

#define HECMW_PARTLEX_V_ND 3015

#define HECMW_PARTLEX_V_USER 3016

#define HECMW_PARTLEX_V_DEFAULT 3052

#define HECMW_PARTLEX_V_AGGREGATE 3053

#define HECMW_PARTLEX_V_DISTRIBUTE 3054

#define HECMW_PARTLEX_V_SIMPLE 3055

extern double HECMW_partlex_get_number(void);

extern char *HECMW_partlex_get_text(void);

extern int HECMW_partlex_get_lineno(void);

extern int HECMW_partlex_next_token(void);

extern int HECMW_partlex_next_token_skip(int skip_token);

extern int HECMW_partlex_set_input(FILE *fp);

extern int HECMW_partlex_skip_line(void);

extern int HECMW_partlex_unput_token(void);

#endif
