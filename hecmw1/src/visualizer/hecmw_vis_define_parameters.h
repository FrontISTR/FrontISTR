#ifndef HECMW_VIS_DEFINE_PARAMETERS_H_INCLUDED
#define HECMW_VIS_DEFINE_PARAMETERS_H_INCLUDED

void transform_range_vertex(double range[6], double vertex[24]);
void get_frame_transform_matrix(double view_point_d[3], double screen_point[3], double up[3], double coff_matrix[3][3]);
void find_inverse_matrix(double coff_matrix[3][3], double inv_matrix[3][3]);
void transform_frame(double screen_point[3], double vertex[24], double coff_matrix[3][3],
		double n_vertex[24]);
void transform_frame3(double screen_point[3], double f[3][3], double coff_matrix[3][3],
		double n_f[3][3]);
void transform2_frame(double coff_matrix[3][3], double view_point[3]);
void tranverse_transform(double screen_point[3], double point_s[3], double inv_matrix[3][3], double point_o[3]);
void transform_frame4(double screen_point[3], double iso_p[6], double coff_matrix[3][3], double n_iso[6]);
void find_projection_range3(double view_point[3],double n_iso[6], double pixel_d[2][2], double iso_p[6]);
void find_projection_range2(double view_point[3],  double n_f[3][3],
		double scr_area[4]);
void find_projection_range(double view_point[3],  double n_vertex[24],
		double scr_area[4]);
void view_parameter_define(int ii, int num_of_frames, int rotate_style, double view_point_d[3], double screen_point[3], double up[3],
		int num_of_lights, double *light_point , double trange[6]);
void view1_parameter_define(int ii, int num_of_frames, int rotate_style, double view_point_d[3], double screen_point[3], int num_of_lights, double *light_point, double up[3], double trange[6]);

#endif /* HECMW_VIS_DEFINE_PARAMETERS_H_INCLUDED */
