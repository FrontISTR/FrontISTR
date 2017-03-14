/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_vis_subimage_composite_sf.h"

void composite_subimage_sf(int pesize, int pixn, double *n_subimage, double *n_subopa,
		double *subimage, double *subopa)
{
	int i, j;

	for(j=0;j<pixn;j++) {
		subimage[j*3+0]=n_subimage[j*3+0];
		subimage[j*3+1]=n_subimage[j*3+1];
		subimage[j*3+2]=n_subimage[j*3+2];
		subopa[j]=n_subopa[j];
	}
	for(i=1;i<pesize;i++) {
		for(j=0;j<pixn;j++) {
			if(n_subopa[i*pixn+j]<subopa[j]) {
				subopa[j]=n_subopa[i*pixn+j];
				subimage[j*3]=n_subimage[i*pixn*3+j*3+0];
				subimage[j*3+1]=n_subimage[i*pixn*3+j*3+1];
				subimage[j*3+2]=n_subimage[i*pixn*3+j*3+2];
			}
		}
	}
	return;
}
