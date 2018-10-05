/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
#include "hecmw_fstr_output_vtk.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "hecmw_malloc.h"
#include "hecmw_etype.h"
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_comm_util.h"
#include "hecmw_vis_combine.h"
#include "hecmw_fstr_endian.h"

void vtk_output (struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, char *outfile, char *outfile1, int *max_timestep, HECMW_Comm VIS_COMM)
{
	int i, j, k;
	int jS, jE;
	int myrank, petot, steptot;
	int n_node, n_elem, shift, etype;
	int data_tot_n, data_tot_e;
	int table342[10] = {0, 1, 2, 3, 6, 4, 5, 7, 8, 9};
	char file_pvd[HECMW_FILENAME_LEN], file_pvtu[HECMW_FILENAME_LEN], file_vtu[HECMW_FILENAME_LEN], buf[HECMW_FILENAME_LEN];
	char *data_label;
	static int is_first=0;
	FILE *outfp;
	HECMW_Status stat;

	HECMW_Comm_rank (VIS_COMM, &myrank);
	HECMW_Comm_size (VIS_COMM, &petot);
	n_node = mesh->n_node;
	n_elem = mesh->n_elem;
	data_tot_n = 0;
	for(i=0; i<data->nn_component; i++){
		data_tot_n += data->nn_dof[i];
	}
	data_tot_e = 0;
	for(i=0; i<data->ne_component; i++){
		data_tot_e += data->ne_dof[i];
	}

	sprintf(file_vtu, "%s/%s.%d.vtu", outfile1, outfile, myrank);
	if(HECMW_ctrl_make_subdir(file_vtu)) {
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0009: Cannot open output directory");
	}

	/* outpu pvd file */
	/* if (myrank == 0 && is_first == 0) {
		sprintf(file_pvd, "%s.pvd", outfile1);
		outfp = fopen (file_pvd, "w");
		fprintf (outfp, "<?xml version=\"1.0\"?>\n");
		fprintf (outfp, "<VTKFile type=\"Collection\" version=\"1.0\">\n");
		fprintf (outfp, "<Collection>\n");
		for(i=0; i < *max_timestep+1 ;i++){
			fprintf (outfp, "<DataSet part=\"0\" timestep=\"%d\" file=\"%s.%04d.pvtu\"/>\n", i, outfile, i);
		}
		fprintf (outfp, "</Collection>\n");
		fprintf (outfp, "</VTKFile>\n");
		fclose (outfp);
		is_first = 1;
	}*/

	if (myrank == 0) {
		/* outpu pvtu file */
		sprintf(file_pvtu, "%s.pvtu", outfile1);
		outfp = fopen (file_pvtu, "w");
		fprintf (outfp, "<?xml version=\"1.0\"?>\n");
		fprintf (outfp, "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"%s\">\n", HECMW_endian_str());
		fprintf (outfp, "<PUnstructuredGrid>\n");
		fprintf (outfp, "<PPoints>\n");
		fprintf (outfp, "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
		fprintf (outfp, "</PPoints>\n");
		fprintf (outfp, "<PCells>\n");
		fprintf (outfp, "<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>\n");
		fprintf (outfp, "<PDataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>\n");
		fprintf (outfp, "<PDataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"/>\n");
		fprintf (outfp, "</PCells>\n");
		fprintf (outfp, "<PPointData>\n");
		for(i=0; i<data->nn_component; i++){
			fprintf (outfp, "<PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"/>\n", data->node_label[i], data->nn_dof[i]);
		}
		fprintf (outfp, "</PPointData>\n");
		fprintf (outfp, "<PCellData>\n");
		fprintf (outfp, "<PDataArray type=\"Int16\" Name=\"Mesh_Type\" NumberOfComponents=\"1\" format=\"ascii\"/>\n");
		for(i=0; i<data->ne_component; i++){
			fprintf (outfp, "<PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"/>\n", data->elem_label[i], data->ne_dof[i]);
		}
		fprintf (outfp, "</PCellData>\n");
		for(i=0; i<petot; i++){
			sprintf (buf, "./%s/%s.%d.vtu", outfile, outfile, i);
			fprintf (outfp, "<Piece Source=\"%s\"/>\n", buf);
		}
		fprintf (outfp, "</PUnstructuredGrid>\n");
		fprintf (outfp, "</VTKFile>\n");
		fclose (outfp);
	}

	/* outpu vtu file */
	outfp = fopen (file_vtu, "w");
	fprintf (outfp, "<?xml version=\"1.0\"?>\n");
	fprintf (outfp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n");
	fprintf (outfp, "<UnstructuredGrid>\n");
	fprintf (outfp, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_node, n_elem);
	fprintf (outfp, "<Points>\n");
	fprintf (outfp, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for(i=0; i<n_node; i++){
		fprintf (outfp, "%e %e %e\n", (float)mesh->node[3*i], (float)mesh->node[3*i+1], (float)mesh->node[3*i+2]);
	}
	fprintf (outfp, "</DataArray>\n");
	fprintf (outfp, "</Points>\n");
	fprintf (outfp, "<Cells>\n");
	fprintf (outfp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
	for(i=0; i<n_elem; i++){
		jS=mesh->elem_node_index[i];
		jE=mesh->elem_node_index[i+1];
		shift=0;
		if(mesh->elem_type[i]==641) shift=2;
		if(mesh->elem_type[i]==761) shift=3;
		if(mesh->elem_type[i]==781) shift=4;
		if(mesh->elem_type[i]==342){
			for(j=jS; j<jE-shift; j++){
				fprintf (outfp, "%d ", mesh->elem_node_item[jS+table342[j-jS]]-1);
			}
		}else{
			for(j=jS; j<jE-shift; j++){
				fprintf (outfp, "%d ", mesh->elem_node_item[j]-1);
			}
		}
		fprintf (outfp, "\n");
	}
	fprintf (outfp, "</DataArray>\n");
	fprintf (outfp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
	shift=0;
	for(i=0; i<n_elem; i++){
		if(mesh->elem_type[i]==641) shift+=2;
		if(mesh->elem_type[i]==761) shift+=3;
		if(mesh->elem_type[i]==781) shift+=4;
		fprintf (outfp, "%d ", mesh->elem_node_index[i+1]-shift);
	}
	fprintf (outfp, "\n");
	fprintf (outfp, "</DataArray>\n");
	fprintf (outfp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
	for(i=0; i<n_elem; i++){
		fprintf (outfp, "%d ", HECMW_get_etype_vtk_shape(mesh->elem_type[i]));
	}
	fprintf (outfp, "\n");
	fprintf (outfp, "</DataArray>\n");
	fprintf (outfp, "</Cells>\n");
	fprintf (outfp, "<PointData>\n");
	for(i=0; i<data->nn_component; i++){
		fprintf (outfp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n", data->node_label[i], data->nn_dof[i]);
		shift=0;
		for(j=0; j<i; j++){
			shift += data->nn_dof[j];
		}
		for(j=0; j<n_node; j++){
			for(k=0; k<data->nn_dof[i]; k++){
				fprintf (outfp, "%e ", (float)data->node_val_item[j*data_tot_n+k+shift]);
			}
			fprintf (outfp, "\n");
		}
		fprintf (outfp, "</DataArray>\n");
	}
	fprintf (outfp, "</PointData>\n");
	fprintf (outfp, "<CellData>\n");
	fprintf (outfp, "<DataArray type=\"Int16\" Name=\"Mesh_Type\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for(i=0; i<n_elem; i++){
		fprintf (outfp, "%d ", mesh->elem_type[i]);
	}
	fprintf (outfp, "\n");
	fprintf (outfp, "</DataArray>\n");
	for(i=0; i<data->ne_component; i++){
		fprintf (outfp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n", data->elem_label[i], data->ne_dof[i]);
		shift=0;
		for(j=0; j<i; j++){
			shift += data->ne_dof[j];
		}
		for(j=0; j<n_elem; j++){
			for(k=0; k<data->ne_dof[i]; k++){
				fprintf (outfp, "%e ", (float)data->elem_val_item[j*data_tot_e+k+shift]);
			}
			fprintf (outfp, "\n");
		}
		fprintf (outfp, "</DataArray>\n");
	}
	fprintf (outfp, "</CellData>\n");
	fprintf (outfp, "</Piece>\n");
	fprintf (outfp, "</UnstructuredGrid>\n");
	fprintf (outfp, "</VTKFile>\n");
	fclose (outfp);
}

void bin_vtk_output (struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, char *outfile, char *outfile1, int *max_timestep, HECMW_Comm VIS_COMM)
{
	int i, j, k;
	int jS, jE;
	int myrank, petot, steptot;
	int n_node, n_elem, shift, etype;
	int data_tot_n, data_tot_e, in, ioffset;
	int *offset;
	uint8_t uint8;
	uint16_t uint16;
	uint32_t uint32;
	uint64_t uint64;
	float val, val1, val2, val3;
	char file_pvd[HECMW_FILENAME_LEN], file_pvtu[HECMW_FILENAME_LEN], file_vtu[HECMW_FILENAME_LEN], buf[HECMW_FILENAME_LEN];
	char *data_label;
	static int is_first=0;
	int table342[10] = {0, 1, 2, 3, 6, 4, 5, 7, 8, 9};
	FILE *outfp;
	HECMW_Status stat;

	HECMW_Comm_rank (VIS_COMM, &myrank);
	HECMW_Comm_size (VIS_COMM, &petot);
	n_node = mesh->n_node;
	n_elem = mesh->n_elem;
	data_tot_n = 0;
	for(i=0; i<data->nn_component; i++){
		data_tot_n += data->nn_dof[i];
	}
	data_tot_e = 0;
	for(i=0; i<data->ne_component; i++){
		data_tot_e += data->ne_dof[i];
	}

	sprintf(file_vtu,  "%s/%s.%d.vtu", outfile1, outfile, myrank);
	if(HECMW_ctrl_make_subdir(file_vtu)) {
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0009: Cannot open output directory");
	}

	/* output pvd file */
	/* if (myrank == 0 && is_first == 0) {
		sprintf(file_pvd,  "%s.pvd",  outfile1);
		outfp = fopen (file_pvd, "w");
		fprintf (outfp, "<?xml version=\"1.0\"?>\n");
		fprintf (outfp, "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"%s\" header_type=\"UInt32\">\n", HECMW_endian_str());
		fprintf (outfp, "<Collection>\n");
		for(i=0; i<*max_timestep ;i++){
			fprintf (outfp, "<DataSet part=\"0\" timestep=\"%d\" file=\"mesh_vis_psf.%04d.pvtu\"/>\n", i+1, i+1);
		}
		fprintf (outfp, "</Collection>\n");
		fprintf (outfp, "</VTKFile>\n");
		fclose (outfp);
	} */

	if (myrank == 0) {
		/* outpu pvtu file */
		sprintf(file_pvtu, "%s.pvtu", outfile1);
		outfp = fopen (file_pvtu, "w");
		fprintf (outfp, "<?xml version=\"1.0\"?>\n");
		fprintf (outfp, "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"%s\" header_type=\"UInt32\">\n", HECMW_endian_str());
		fprintf (outfp, "<PUnstructuredGrid>\n");
		fprintf (outfp, "<PPoints>\n");
		fprintf (outfp, "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
		fprintf (outfp, "</PPoints>\n");
		fprintf (outfp, "<PCells>\n");
		fprintf (outfp, "<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\"/>\n");
		fprintf (outfp, "<PDataArray type=\"Int32\" Name=\"offsets\" format=\"appended\"/>\n");
		fprintf (outfp, "<PDataArray type=\"Int32\" Name=\"types\" format=\"appended\"/>\n");
		fprintf (outfp, "</PCells>\n");
		fprintf (outfp, "<PPointData>\n");
		for(i=0; i<data->nn_component; i++){
			fprintf (outfp, "<PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"appended\"/>\n", data->node_label[i], data->nn_dof[i]);
		}
		fprintf (outfp, "</PPointData>\n");
		fprintf (outfp, "<PCellData>\n");
		for(i=0; i<data->ne_component; i++){
			fprintf (outfp, "<PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"appended\"/>\n", data->elem_label[i], data->ne_dof[i]);
		}
		fprintf (outfp, "<PDataArray type=\"Int32\" Name=\"Mesh_Type\" NumberOfComponents=\"1\" format=\"appended\"/>\n");
		fprintf (outfp, "</PCellData>\n");
		for(i=0; i<petot; i++){
			sprintf (buf,  "./%s/%s.%d.vtu", outfile, outfile, i);
			fprintf (outfp, "<Piece Source=\"%s\"/>\n", buf);
		}
		fprintf (outfp, "</PUnstructuredGrid>\n");
		fprintf (outfp, "</VTKFile>\n");
		fclose (outfp);
		is_first = 1;
	}

	/* outpu vtu file */
	ioffset = 5 + data->nn_component + data->ne_component;
	offset = HECMW_malloc(sizeof(int)*ioffset);

	uint64 = 0;
	for(i=0; i<n_elem; i++){
		jS=mesh->elem_node_index[i];
		jE=mesh->elem_node_index[i+1];
		shift=0;
		if(mesh->elem_type[i]==641) shift=2;
		if(mesh->elem_type[i]==761) shift=3;
		if(mesh->elem_type[i]==781) shift=4;
		uint64 += jE-shift-jS;
	}

	offset[0] = 0;
	offset[1] = offset[0] + sizeof(int) + 3*n_node   *sizeof(float);
	offset[2] = offset[1] + sizeof(int) + (int)uint64*sizeof(int);
	offset[3] = offset[2] + sizeof(int) + n_elem     *sizeof(int);
	offset[4] = offset[3] + sizeof(int) + n_elem     *sizeof(int);
	for(i=0; i<data->nn_component; i++){
		offset[5+i] = offset[4+i] + sizeof(int) + data->nn_dof[i]*n_node*sizeof(int);
	}
	for(i=0; i<data->ne_component; i++){
		offset[5+data->nn_component+i] = offset[4+data->nn_component+i] + sizeof(int) + data->ne_dof[i]*n_elem*sizeof(int);
	}

	outfp = fopen (file_vtu, "w");
	fprintf (outfp, "<?xml version=\"1.0\"?>\n");
	fprintf (outfp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"%s\" header_type=\"UInt32\">\n", HECMW_endian_str());
	fprintf (outfp, "<UnstructuredGrid>\n");
	fprintf (outfp, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_node, n_elem);
	fprintf (outfp, "<Points>\n");
	fprintf (outfp, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\">\n", offset[0]);
	fprintf (outfp, "</DataArray>\n");
	fprintf (outfp, "</Points>\n");
	fprintf (outfp, "<Cells>\n");
	fprintf (outfp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%d\">\n", offset[1]);
	fprintf (outfp, "</DataArray>\n");
	fprintf (outfp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%d\">\n", offset[2]);
	fprintf (outfp, "</DataArray>\n");
	fprintf (outfp, "<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%d\">\n", offset[3]);
	fprintf (outfp, "</DataArray>\n");
	fprintf (outfp, "</Cells>\n");
	fprintf (outfp, "<PointData>\n");

	for(i=0; i<data->nn_component; i++){
		fprintf (outfp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"appended\" offset=\"%d\">\n", data->node_label[i], data->nn_dof[i], offset[4+i]);
		fprintf (outfp, "</DataArray>\n");
	}

	fprintf (outfp, "</PointData>\n");
	fprintf (outfp, "<CellData>\n");
	for(i=0; i<data->ne_component; i++){
		fprintf (outfp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"appended\" offset=\"%d\">\n", data->elem_label[i], data->ne_dof[i], offset[4+data->nn_component+i]);
		fprintf (outfp, "</DataArray>\n");
	}

	fprintf (outfp, "<DataArray type=\"Int32\" Name=\"Mesh_Type\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n", offset[ioffset-1]);
	fprintf (outfp, "</DataArray>\n");
	fprintf (outfp, "</CellData>\n");
	fprintf (outfp, "</Piece>\n");
	fprintf (outfp, "</UnstructuredGrid>\n");
	fprintf (outfp, "<AppendedData encoding=\"raw\">\n");

	fprintf (outfp, " _");
	uint32 = (uint32_t)(3*n_node*sizeof(float));
	fwrite (&uint32, sizeof(uint32), 1, outfp);
	for(i=0; i<n_node; i++){
		val = (float)mesh->node[3*i];
		fwrite (&val, sizeof(float), 1, outfp);
		val = (float)mesh->node[3*i+1];
		fwrite (&val, sizeof(float), 1, outfp);
		val = (float)mesh->node[3*i+2];
		fwrite (&val, sizeof(float), 1, outfp);
	}

	uint32 = (uint32_t)(uint64*sizeof(int));
	fwrite (&uint32, sizeof(uint32), 1, outfp);
	for(i=0; i<n_elem; i++){
		jS=mesh->elem_node_index[i];
		jE=mesh->elem_node_index[i+1];
		shift=0;
		if(mesh->elem_type[i]==641) shift=2;
		if(mesh->elem_type[i]==761) shift=3;
		if(mesh->elem_type[i]==781) shift=4;
		if(mesh->elem_type[i]==342){
			for(j=jS; j<jE-shift; j++){
				in = (int)mesh->elem_node_item[jS+table342[j-jS]]-1;
				fwrite (&in, sizeof(int), 1, outfp);
			}
		}else{
			for(j=jS; j<jE-shift; j++){
				in = (int)mesh->elem_node_item[j]-1;
				fwrite (&in, sizeof(int), 1, outfp);
			}
		}
	}

	uint32 = (uint32_t)(n_elem*sizeof(int));
	fwrite (&uint32, sizeof(uint32), 1, outfp);
	shift=0;
	for(i=0; i<n_elem; i++){
		if(mesh->elem_type[i]==641) shift+=2;
		if(mesh->elem_type[i]==761) shift+=3;
		if(mesh->elem_type[i]==781) shift+=4;
		in = (int)mesh->elem_node_index[i+1]-shift;
		fwrite (&in, sizeof(int), 1, outfp);
	}

	uint32 = (uint32_t)(n_elem*sizeof(int));
	fwrite (&uint32, sizeof(uint32), 1, outfp);
	for(i=0; i<n_elem; i++){
		//uint8 = (uint8_t)HECMW_get_etype_vtk_shape(mesh->elem_type[i]);
	  //fwrite (&uint8, sizeof(u_int8_t), 1, outfp);
		in = (int)HECMW_get_etype_vtk_shape(mesh->elem_type[i]);
	  fwrite (&in, sizeof(int), 1, outfp);
	}

	for(i=0; i<data->nn_component; i++){
		uint32 = (uint32_t)(data->nn_dof[i]*n_node*sizeof(int));
		fwrite (&uint32, sizeof(uint32), 1, outfp);

		shift=0;
		for(j=0; j<i; j++){
			shift += data->nn_dof[j];
		}
		for(j=0; j<n_node; j++){
			for(k=0; k<data->nn_dof[i]; k++){
				val = (float)data->node_val_item[j*data_tot_n+k+shift];
				fwrite (&val, sizeof(float), 1, outfp);
			}
		}
	}

	for(i=0; i<data->ne_component; i++){
		uint32 = (uint32_t)(data->ne_dof[i]*n_elem*sizeof(int));
		fwrite (&uint32, sizeof(uint32), 1, outfp);

		shift=0;
		for(j=0; j<i; j++){
			shift += data->ne_dof[j];
		}
		for(j=0; j<n_elem; j++){
			for(k=0; k<data->ne_dof[i]; k++){
				val = (float)data->elem_val_item[j*data_tot_e+k+shift];
				fwrite (&val, sizeof(float), 1, outfp);
			}
		}
	}

	uint32 = (uint32_t)(n_elem*sizeof(int));
	fwrite (&uint32, sizeof(uint32), 1, outfp);
	for(i=0; i<n_elem; i++){
		//uint16 = (uint16_t)mesh->elem_type[i];
	  //fwrite (&uint16, sizeof(u_int16_t), 1, outfp);
		in = (int)mesh->elem_type[i];
	  fwrite (&in, sizeof(int), 1, outfp);
	}

	fprintf (outfp, "</AppendedData>\n");
	fprintf (outfp, "</VTKFile>\n");
	fclose (outfp);
}

void HECMW_vtk_output (struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, char *outfile, char *outfile1, int *max_timestep, HECMW_Comm VIS_COMM)
{
	vtk_output (mesh, data, outfile, outfile1, max_timestep, VIS_COMM);
}

void HECMW_bin_vtk_output (struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, char *outfile, char *outfile1, int *max_timestep, HECMW_Comm VIS_COMM)
{
	bin_vtk_output (mesh, data, outfile, outfile1, max_timestep, VIS_COMM);
}
