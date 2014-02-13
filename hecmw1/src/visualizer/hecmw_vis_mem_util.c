/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Visualization                                     *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#include "hecmw_vis_mem_util.h"

#include <stdio.h>
#include <stdlib.h>
#include "hecmw_util.h"


void HECMW_vis_memory_exit(char *var)
{
	fprintf(stderr, "#### HEC-MW-VIS-E0001:There is no enough memory allocated for variable %s\n", var);
	HECMW_finalize();
	exit(0);
}

void HECMW_vis_print_exit(char *var)
{
	fprintf(stderr, "%s\n", var);
	HECMW_finalize();
	exit(0);
}


void mfree(void *pointer)
{
	HECMW_free(pointer);
	pointer = NULL;
}

Point *alloc_verts(int num)
{
	int i;
	Point *verts;

	if ((verts = (Point *) HECMW_calloc(num,sizeof(Point))) == NULL) {
		fprintf(stderr,"There is not enough memory, alloc_verts\n");
		return NULL;
	}

	for (i = 0; i < (num - 1); i++) {
		(verts + i)->nextpoint = (verts + i + 1);
		(verts + i)->ident = 0;
	}
	(verts + num - 1)->ident = 0;
	(verts + num - 1)->nextpoint = NULL;


	return verts;
}


Polygon *alloc_polygons(int num)
{
	int	i;
	Polygon *polygons;

	if ((polygons = (Polygon *) HECMW_calloc(num,sizeof(Polygon))) == NULL) {
		fprintf(stderr,"There is not enough memory, alloc_polygons\n");
		return NULL;
	}

	for (i = 0; i < (num - 1); i++) {
		(polygons + i)->nextpolygon = (polygons + i + 1);
		(polygons + i)->plist = NULL;
	}
	(polygons + num - 1)->nextpolygon = NULL;
	(polygons + num - 1)->plist = NULL;


	return polygons;
}
