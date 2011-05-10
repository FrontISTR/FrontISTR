#ifndef HECMW_VIS_FONT_TEXTURE_H_INCLUDED
#define HECMW_VIS_FONT_TEXTURE_H_INCLUDED

void font7_generate(char input, int output[7][7]);
void mark_time_label(double font_size,
		int xr, int yr,
		double font_color[3], double background_color[3],
		double start_time, double time_interval,
		int timestep, int max_len_step, double *image);
void value2_to_rgb(double value, double color[3], int color_system_type);
void generate_color_bar(int scale_marking_on, double font_size, int color_bar_style,
		int mark_0_on, int color_mapping_bar_on,
		int xr, int yr, double font_color[3], int color_system_type,
		int color_mapping_style, double *interval_point,
		int interval_mapping_num, int num_of_scale,
		double tmincolor, double tmaxcolor,
		double org_mincolor, double org_maxcolor, double *image);
unsigned short int change_short_int_order(unsigned short int n);
unsigned int change_unsigned_int_order(unsigned int n);
int change_int_order(int n);

#endif /* HECMW_VIS_FONT_TEXTURE_H_INCLUDED */
