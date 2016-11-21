/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_LIB_FC_INCLUDED
#define HECMW_LIB_FC_INCLUDED


extern char *HECMW_strcpy_f2c(const char *fstr, int flen);


extern char *HECMW_strcpy_f2c_r(const char *fstr, int flen, char *buf, int bufsize);


extern int HECMW_strcpy_c2f(const char *cstr, char *fstr, int flen);

#endif
