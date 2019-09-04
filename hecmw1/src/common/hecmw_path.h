/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_PATH_INCLUDED
#define HECMW_PATH_INCLUDED

extern int HECMW_get_path_separator(void);

extern int HECMW_is_absolute_path(const char *path);

extern char *HECMW_dirname(const char *path);

extern char *HECMW_basename(const char *path);

#endif
