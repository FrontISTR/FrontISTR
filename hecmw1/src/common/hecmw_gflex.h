/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_GFLEX_INCLUDED
#define HECMW_GFLEX_INCLUDED

#include <stdio.h>


enum {
	HECMW_GFLEX_NL = 1000,
	HECMW_GFLEX_INT,
	HECMW_GFLEX_DOUBLE,
	HECMW_GFLEX_NAME,
};


int HECMW_gflex_get_lineno(void);


double HECMW_gflex_get_number(void);


char *HECMW_gflex_get_text(void);


int HECMW_gflex_next_token(void);


int HECMW_gflex_next_token_skip(int skip_token);


int HECMW_gflex_set_input(FILE *fp);


int HECMW_gflex_skip_line(void);

#endif
