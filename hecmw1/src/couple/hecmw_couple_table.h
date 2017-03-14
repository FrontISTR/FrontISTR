/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_TABLE
#define INC_HECMW_COUPLE_TABLE


static int hecmw_surf_node_table_tet1[4][3] = {
    { 2, 3, 4 },
    { 1, 4, 3 },
    { 1, 2, 4 },
    { 1, 3, 2 },
};


static int hecmw_surf_node_table_pri1[5][4] = {
    { 2, 3, 6, 5 },
    { 3, 1, 4, 6 },
    { 1, 2, 5, 4 },
    { 3, 2, 1 },
    { 4, 5, 6 },
};


static int hecmw_surf_node_table_hex1[6][4] = {
    { 4, 1, 5, 8 },
    { 2, 3, 7, 6 },
    { 1, 2, 6, 5 },
    { 3, 4, 8, 7 },
    { 4, 3, 2, 1 },
    { 5, 6, 7, 8 },
};

#endif	/* INC_HECMW_COUPLE_TABLE */
