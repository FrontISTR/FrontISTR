/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_CONN_CONV_INCLUDED
#define HECMW_CONN_CONV_INCLUDED

#define HECMW_CONNTYPE_HECMW 1
#define HECMW_CONNTYPE_ABAQUS 2
#define HECMW_CONNTYPE_NASTRAN 3


extern int HECMW_convert_connectivity(int from, int hecmw_etype, int *conn);

#endif
