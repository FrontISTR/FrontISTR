#ifndef HECMW_VIS_READ_CONTROL_H_INCLUDED
#define HECMW_VIS_READ_CONTROL_H_INCLUDED

#include <stdio.h>
#include "hecmw_vis_SF_geom.h"
#include "hecmw_vis_ray_trace.h"

#define MAX_LINE_LEN  256

int is_blank_line(char *buf);
int is_comment_line(char *buf);
void get_string_item(char *para, char *buf, int *start_location, char para2[128]);
int get_int_item(char *para, char *buf, int *start_location);
double get_double_item(char *para, char *buf, int *start_location);
int get_keyword_item(char *buf, char *para);
void HECMW_vis_read_control(FILE *fp, int pesize, int mynode, PSF_link *psf, PVR_link *pvr);

#endif /* HECMW_VIS_READ_CONTROL_H_INCLUDED */


