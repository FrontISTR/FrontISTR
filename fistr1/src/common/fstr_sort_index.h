/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/**
 * @brief I/O and Utility
 */

/*
        Index Sorting for FSTR
        2004.10.18 by N.Imai
        -------------------------
        [Fortran]
        integer(kind=4) :: index_data(2,:), n
        call fstr_sort_index( index_data, n )
*/

#ifndef fstr_sort_indexH
#define fstr_sort_indexH

void c_fstr_sort_index(int* index_data, int n);

#endif
