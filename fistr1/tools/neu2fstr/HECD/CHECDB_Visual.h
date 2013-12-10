/*=====================================================================*
 *                                                                     *
 *   Software Name : neu2fstr                                          *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : FrontSTR Input File Converter                     *
 *                                                                     *
 *            Written by Noboru Imai (Univ. of Tokyo)                  *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

/*
	CHECDB_Visual Ver.1.0
*/


#ifndef CHECDB_VisualH
#define CHECDB_VisualH

#include "CHECDataBlock.h"


// parameters for screen
class CVis_ViewParams {
public:
	// screen

	int x_resolution;		// default:512
	int y_resolution;		// default:512
	int num_of_lights;		// default:1
	double* position_of_lights;	// (default:NULL--not specified)
	double viewpoint[3];
	bool fg_viewpoint;		// (inner use, default:false);
	double look_at_point[3];
	bool fg_look_at_point;		// (inner use, default:false);
	double up_direction[3];		// default:{0,0,1}
	double ambient_coef;		// default:0.3
	double diffuse_coef;		// default:0.7
	double specular_coef;		// default:0.6
	int color_mapping_style;
	enum {
		color_mapping_style_linear = 1, // default;
		color_mapping_style_clipped_linear,
		color_mapping_style_nonlinear,
		color_mapping_style_automatic
	};
	int interval_mapping_num;	// (default:0)
	double* interval_mapping;	// (default:NULL -- not specified)

	int rotate_style;
	enum {
		rotate_style_none = 0, // default ( inner use )
		rotate_style_x,
		rotate_style_y,
		rotate_style_z,
		rotate_style_view_point
	};
	int num_of_frame;		// default:8
	int color_mapping_bar_on;	// default:0
	int scale_marking_on;		// default:0
	int num_of_scales;		// default:3
	double font_size;		// default:1
	double font_color[3];		// default:{1,1,1}
	double background_color[3];	// default:{0,0,0}

	int color_system_type;
	enum {
		color_system_type_blue_red = 1, // default
		color_system_type_rainbow,
		color_system_type_block_white
	};
	int fixed_range_on;		// default:0
	double range_value[2];
	bool fg_range_value;		// (inner use, default:false)

public:
	// methods
	CVis_ViewParams();
	virtual ~CVis_ViewParams();
	virtual void Init();
	virtual void Write( class CHECData* hecd );
	static void WriteVisParam(class CHECData* hecd, const char* name, const char* fmt, ... );
	static void WriteVisPArry(class CHECData* hecd, const char* name, char type, int n, void* p ); // for array

	template <class T>
	static int ReadVisValueT( char* s, int n, T* value, const char* fmt );
	static int ReadVisValue( char* s, int n, int* value );
	static int ReadVisValue( char* s, int n, double* value );
	static int ReadVisValue( char* s, int n, char* value ); // n:sizeof(value), return:strlen(value)
	template <class T>
	static T* ReadVisPArry( char* s, int& n, const char* fmt );

	virtual bool ReadData( const char* line );
	virtual bool Read( class CHECData* hecd );
};


//-----------------------------------------------------------------------------

class CVis_PSR  : public CVis_ViewParams {
public:
	// common ------------------------------------------------

	int surface_style;
	enum {
		surface_style_boundary = 1, // default
		surface_style_equivalent,
		surface_style_user
	};
	int display_method;
	enum {
		display_method_color = 1, // default
		display_method_boundary,
		display_method_color_and_boundary,
		display_method_mono_tone,
		display_method_contour
	};
	char color_comp_name[100];
	char color_subcomp_name[4];	// "norm", "x"(default), "y" or "z"
	int color_comp;			// default:0
	int color_subcomp;		// default:1
	int iso_number;			// default:5
	double specified_color; 	// effective if display_method == 4
	int deform_display_on;		// 1:on, 0:off, default:0
	char deform_comp_name[100];
	int deform_comp;		// default:0
	double deform_scale;
	bool fg_deform_scale; // is deform_scale set? (inner use) default:false
	int initial_style;
	enum {
		initial_style_none = 0, // default
		initial_style_mesh,
		initial_style_fill,
		initial_style_shading,
		initial_style_dot_mesh
	};
	int deform_style;
	enum {
		deform_style_none = 0,
		deform_style_mesh,
		deform_style_fill,
		deform_style_shading,
		deform_style_dot_mesh // default
	};
	double initial_line_color[3];	// default:{0,0,1}
	double deform_line_color[3];	// default:{0,0,1}
	char output_type[3];	// "AVS"(default), "BMP"

	// for surface_style == surface_style_equivalent ---------

	char data_comp_name[100];
	char data_subcomp_name[4];	// "norm", "x"(default), "y" or "z"
	int data_comp;			// default:0
	int data_subcomp;		// default:1
	double iso_value;

	// for surface_style == surface_style_user ---------------

	int method;
	enum {
		method_sphere = 1,
		method_ellipsoid,
		method_hyperbola,
		method_parabola,
		method_quadric		// default
	};
	double point[3];		// default:{0,0,0}
	double radius;			// default:1.0
	double length;			// (default:1.0)
	double coef[10];		// (default:0.0)

	// rendering parameters(output_type == "BMP") ------------
	// parameters defined in CVis_ScreenParams and followings

	double isoline_color[3];
	bool fg_isoline_color;		// (inner use, default:false )
	int boundary_line_on;		// default:0

public:
	// methods
	CVis_PSR();
	virtual ~CVis_PSR();
	virtual void Init();
	virtual void Write( class CHECData* hecd );
	virtual bool ReadData( const char* line );
};


//-----------------------------------------------------------------------------

class CVis_PVR : public CVis_ViewParams {
public:
	// for dividing ------------------------------------------

	int maximum_refinement;		// default:100
	int n_voxel_x;
	bool fg_n_voxel_x;		// (inner use, default:false)
	int n_voxel_y;
	bool fg_n_voxel_y;		// (inner use, default:false)
	int n_voxel_z;
	bool fg_n_voxel_z;		// (inner use, default:false)
	char voxel_filename[100];	// default:""
	int x_specified_level;		// default:100
	int y_specified_level;		// default:100
	int z_specified_level;		// default:100

	// for opacity

	int transfer_function_type;	// default:1
	double opa_value;		// default:0.02
	int num_of_features;		// (default:0)
	double* fea_point;		// (default:NULL (no set))
	char name_lookup[100];		// (default:"" (no set))

	// parameters defined in CVis_ViewParams and followings

	int histogram_on;		// default:0
	double display_range[6];	// (default:0)
	bool fg_display_range;		// (inner use, default:false)

public:
	// methods
	CVis_PVR();
	virtual ~CVis_PVR();
	virtual void Init();
	virtual void Write( class CHECData* hecd );
	virtual bool ReadData( const char* line );
};


//-----------------------------------------------------------------------------
// CHECDB_Visual
//-----------------------------------------------------------------------------


class CHECDB_Visual : public CHECDataBlock {
public:
	int visual_start_step;
	int visual_interval_step;
	int visual_end_step;

 	int surface_num;
	std::vector<CVis_PSR*> psr; // psr.size() == surface_num

	CVis_PVR* pvr; // single instance

	CHECDB_Visual();
	virtual ~CHECDB_Visual();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


#endif

