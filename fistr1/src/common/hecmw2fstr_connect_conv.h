/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

/**
 * HECMW to FSTR Mesh Data Converter
 * Convering Conectivity of Element Type 232, 342 and 352 
 */

void c_hecmw2fstr_connect_conv(int n_elem, int elem_type[],
                               int elem_node_index[], int elem_node_item[]);

/**
 * FSTR to HECMW Mesh Data Converter
 * Convering Conectivity of Element Type 232, 342 and 352
 */
void c_fstr2hecmw_connect_conv(int n_elem, int elem_type[],
                               int elem_node_index[], int elem_node_item[]);
void c_fstr2hecmw_elem_conv(int elem_type, int node[]);
