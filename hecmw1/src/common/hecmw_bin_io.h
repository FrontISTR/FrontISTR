/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef hecmw_bin_ioH
#define hecmw_bin_ioH

#include <stdio.h>

/*---------------------------------------------------------------------------*/
void bfwrite(const void *data, size_t size, size_t nitems, FILE* stream);
void bfwrite_flush(FILE* stream);

/*---------------------------------------------------------------------------*/
/* CAUTION) hecmw_set_endian_info must be executed before calling following
 * functions. */

void hecmw_set_endian_info(void);

/*---------------------------------------------------------------------------*/

int hecmw_write_bin_value(unsigned char* x, int size, FILE* fp);
int hecmw_write_bin_int(int x, FILE* fp);
int hecmw_write_bin_int_arr(int* x, int n, FILE* fp);
int hecmw_write_bin_double(double x, FILE* fp);
int hecmw_write_bin_double_arr(double* x, int n, FILE* fp);

/*---------------------------------------------------------------------------*/

int hecmw_read_bin_value(unsigned char* x, int size, FILE* fp);
int hecmw_read_bin_int(int* x, FILE* fp);
int hecmw_read_bin_int_arr(int* x, int n, FILE* fp);
int hecmw_read_bin_double(double* x, FILE* fp);
int hecmw_read_bin_double_arr(double* x, int n, FILE* fp);

/*---------------------------------------------------------------------------*/

/*
 * fmt : type and array size of arguments
 *       Each cherecter in fmt specifies type of argument and
 *       number after the cheracter means array size.
 * meening of character in fmt
 * 'I' : int
 * 'F' : double
 * 'S' : string (char*)
 * exp)
 *     int n, i[10]; char s[20] = "123";
 *     hecmw_write_bin(fp, "II10S", n, i, s);
 */

int hecmw_write_bin(FILE* fp, const char* fmt, ...);
int hecmw_read_bin(FILE* fp, const char* fmt, ...);

#endif
