/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_ERROR_INCLUDED
#define HECMW_ERROR_INCLUDED

#include <stdarg.h>

extern int HECMW_set_verror(int errorno, const char *fmt, va_list ap);

extern int HECMW_set_error(int errorno, const char *fmt, ...);

extern int HECMW_get_error(char **errmsg);

extern int HECMW_get_errno(void);

extern char *HECMW_get_errmsg(void);

#endif
