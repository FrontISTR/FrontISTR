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


#include <stdarg.h>
#include <vector>
#include "CHECDB.h"
#include "CHECDB_Visual.h"
#include "CHECData.h"
#include "hecd_util.h"


using namespace std;
using namespace hecd_util;


//*****************************************************************************
// CVis_ViewParams
//*****************************************************************************


CVis_ViewParams::CVis_ViewParams()
 : position_of_lights(0), interval_mapping(0)
{
	Init();
}


CVis_ViewParams::~CVis_ViewParams()
{
	delete position_of_lights;
	delete interval_mapping;
}


void CVis_ViewParams::Init()
{
	x_resolution = 512;
	y_resolution = 512;
	num_of_lights = 1;
	delete position_of_lights;
	position_of_lights = 0;
	viewpoint[0] = viewpoint[1] = viewpoint[2] = 0;
	fg_viewpoint = false;
	look_at_point[0] = look_at_point[1] = look_at_point[2] = 0;
	fg_look_at_point = false;
	up_direction[0] = 0;
	up_direction[1] = 0;
	up_direction[2] = 1;
	ambient_coef = 0.3;
	diffuse_coef = 0.7;
	specular_coef = 0.6;
	color_mapping_style = color_mapping_style_linear;
	interval_mapping_num = 0;

	delete[] interval_mapping;
	interval_mapping = 0;

	interval_mapping = 0;
	rotate_style =  rotate_style_none;
	num_of_frame = 8;
	color_mapping_bar_on = 0;
	scale_marking_on = 0;
	num_of_scales = 3;
	font_size = 1;
	font_color[0] = font_color[1] = font_color[2] = 1;
	background_color[0] = background_color[1] = background_color[2] = 0;
	color_system_type = color_system_type_blue_red;
	fixed_range_on = 0;
	range_value[0] = range_value[1] = 0;
	fg_range_value = false;
}



void CVis_ViewParams::WriteVisParam(class CHECData* hecd, const char* name, const char* fmt, ... )
{
	va_list va;
	va_start( va, fmt );
	char buff[256];
	char ss[256];

	sprintf( buff, "!%s", name );

	int n = strlen(fmt);
	if(n>1) {
		strcat(buff, " =");
	}
	for(int i=0; i<n; i++) {
		strcat( buff, " " );
		char c = fmt[i];
		switch(c){
		case 'S':
			strcat( buff, va_arg(va, char*));
			break;
		case 'I':
			sprintf( ss, "%d", va_arg(va, int));
			strcat( buff, ss );
			break;
		case 'F':
			sprintf( ss, "%lg", va_arg(va, double));
			strcat( buff, ss );
			break;
		default:
			assert(0);
		}
	}

	hecd->WriteLine( buff );
	va_end(va);
}


void CVis_ViewParams::WriteVisPArry(class CHECData* hecd, const char* name, char type, int n, void* p )
{
	char buff[256];
	char ss[256];
	int* ip;
	double* dp;

	sprintf( buff, "!%s", name );

	if(n>1) {
		strcat(buff, " =");
	}
	ip = (int*)p;
	dp = (double*)p;
	for(int i=0; i<n; i++) {
		strcat( buff, ", " );
		switch(type){
		case 'I':
			sprintf( ss, "%d", *ip );
			strcat( buff, ss );
			ip++;
			break;
		case 'F':
			sprintf( ss, "%lg", *dp);
			strcat( buff, ss );
			dp++;
			break;
		default:
			assert(0);
		}
	}

	hecd->WriteLine( buff );
}





void CVis_ViewParams::Write( class CHECData* hecd )
{
	WriteVisParam( hecd, "x_resolution",  "I", x_resolution );
	WriteVisParam( hecd, "y_resolution",  "I", y_resolution );
	WriteVisParam( hecd, "num_of_lights", "I", num_of_lights );
	if( num_of_lights > 0 && position_of_lights ){
		WriteVisPArry( hecd, "position_of_lights", 'F', num_of_lights, position_of_lights );
	}
	if( fg_viewpoint )
		WriteVisPArry( hecd, "viewpoint", 'F', 3, viewpoint );
	if( fg_look_at_point )
		WriteVisPArry( hecd, "look_at_point", 'F', 3, look_at_point );
	WriteVisParam( hecd, "up_direction", "F", 3, up_direction );
	WriteVisParam( hecd, "ambient_coef",        "F", ambient_coef );
	WriteVisParam( hecd, "diffuse_coef",        "F", diffuse_coef );
	WriteVisParam( hecd, "specular_coef",       "F", specular_coef );
	WriteVisParam( hecd, "color_mapping_style", "I", color_mapping_style );
	WriteVisParam( hecd, "interval_mapping_num","I", interval_mapping_num );
	if( interval_mapping_num > 0 ) {
		if(!interval_mapping ) assert(0);
		WriteVisPArry( hecd, "nterval_mapping", 'F', interval_mapping_num*2, interval_mapping );
	}
	if( rotate_style != rotate_style_none )
		WriteVisParam( hecd, "rotate_style", "I", rotate_style );
	WriteVisParam( hecd, "num_of_frame", "I", num_of_frame);
	WriteVisParam( hecd, "color_mapping_bar_on", "I", color_mapping_bar_on );
	WriteVisParam( hecd, "scale_marking_on",  "I", scale_marking_on );
	WriteVisParam( hecd, "num_of_scales", "I", num_of_scales );
	WriteVisParam( hecd, "font_size",  "F", font_size );
	WriteVisPArry( hecd, "font_color", 'F', 3, font_color );
	WriteVisPArry( hecd, "font_color", 'F', 3, font_color );
	WriteVisPArry( hecd, "background_color", 'F', 3, background_color );
	WriteVisParam( hecd, "color_system_type", "I", color_system_type );
	WriteVisParam( hecd, "fixed_range_on", "I", fixed_range_on );
	if( fg_range_value )
		WriteVisPArry( hecd, "range_value", 'F', 3, range_value );
}

//-----------------------------------------------------------------------------

template <class T>
int CVis_ViewParams::ReadVisValueT( char* s, int n, T* value, const char* fmt )
{
	const char* delim = ", \t\r\n";
	int i=0;
	char* token = strtok(s, delim);
	while(token && i<n){
		if( sscanf( token, fmt, &value[i]) != 1) return -(i+1);
		token = strtok(0, delim);
		i++;
	}
	return i;
}


int CVis_ViewParams::ReadVisValue( char* s, int n, int* value )
{
	return ReadVisValueT<int>(s, n, value, "%d" );
}


int CVis_ViewParams::ReadVisValue( char* s, int n, double* value )
{
	return ReadVisValueT<double>(s, n, value, "%lf" );
}


int CVis_ViewParams::ReadVisValue( char* s, int n, char* value )
{
	cleanup_token( s );
	int len = strlen(s);
	if( len > n ) return 0;
	strcpy( value, s );
	return len;
}


template <class T>
T* CVis_ViewParams::ReadVisPArry( char* s, int& n, const char* fmt )
{
	const int max_data_n = 200;
	T buff[max_data_n];

	const char* delim = ", \t\r\n";
	n = 0;
	int i=0;
	char* token = strtok(s, delim);
	while(token){
		if( max_data_n >= i) assert(0);
		T x;
		if( sscanf( token, fmt, &x)  != 1) {
			return 0;
		}
		buff[i] = x;
		token = strtok(0, delim);
		i++;
	}

	n = i;
	T* data = new T[n];
	for(i=0; i<n; i++ ) {
		data[i] = buff[i];
	}
	return data;
}


bool CVis_ViewParams::ReadData( const char* line )
{
	char buff[256];
	strcpy( buff, line );
	char* header = strtok( buff, " =,\r\n\t" );
	if( header[0] == 0 || header[0] != '!' ) return false;
	tolower( header );

	#define GENERATE_CODE( x, n, fg ) \
		else if(strcmp( #x, header )==0) {\
			int r = ReadVisValue( 0, n, &x );\
			fg = (r>0); \
			return fg; \
		}
	#define GENERATE_CODEA( x, n, fg ) \
		else if(strcmp( #x, header )==0) {\
			int r = ReadVisValue( 0, n, x );\
			fg = (r>0); \
			return fg; \
		}

	bool dummy;
	if(false); // dummy
	GENERATE_CODE( x_resolution,        1, dummy )
	GENERATE_CODE( y_resolution,        1, dummy )
	GENERATE_CODEA( viewpoint,          3, fg_viewpoint)
	GENERATE_CODEA( look_at_point,      3, fg_look_at_point )
	GENERATE_CODEA( up_direction,       3, dummy )
	GENERATE_CODE( ambient_coef,        1, dummy )
	GENERATE_CODE( diffuse_coef,        1, dummy )
	GENERATE_CODE( specular_coef,       1, dummy )
	GENERATE_CODE( color_mapping_style, 1, dummy )
	GENERATE_CODE( rotate_style,        1, dummy )
	GENERATE_CODE( num_of_frame,        1, dummy )
	GENERATE_CODE( color_mapping_bar_on,1, dummy )
	GENERATE_CODE( scale_marking_on,    1, dummy )
	GENERATE_CODE( num_of_scales,       1, dummy )
	GENERATE_CODE( font_size,           1, dummy )
	GENERATE_CODEA( font_color,         3, dummy )
	GENERATE_CODEA( background_color,   3, dummy )
	GENERATE_CODE( color_system_type,   1, dummy )
	GENERATE_CODE( fixed_range_on,      1, dummy )
	GENERATE_CODEA( range_value,        2, fg_range_value )

	#undef GENERATE_CODE
	#undef GENERATE_CODEA

	#define GENERATE_CODE( x, n, type, fmt ) \
		else if(strcmp( #x, header )==0) {\
			delete[] x; x = 0; \
			x = ReadVisPArry<type>( 0, n, fmt );\
			return x!=0; \
		}

	if(false);
	GENERATE_CODE( position_of_lights, num_of_lights, double, "%lf" )
	GENERATE_CODE( interval_mapping,  interval_mapping_num, double, "%lf" )
	#undef GENERATE_CODE

	return false;
}



bool CVis_ViewParams::Read( CHECData* hecd )
{
	char line[256];
	while(1) {
		if(!hecd->ReadLine( line )) break;
		if(!ReadData( line )) {
			hecd->PushReadLine(line);
			break;
		}
	}
	return true;
}



//*****************************************************************************
// CVis_PSR
//*****************************************************************************



CVis_PSR::CVis_PSR()
 : CVis_ViewParams()
{
	Init();
}


CVis_PSR::~CVis_PSR()
{
}


void CVis_PSR::Init()
{
	// common ------------------------------------------------

	surface_style = surface_style_boundary;
	display_method = display_method_color;
	color_comp_name[0] = 0;
	strcpy( color_subcomp_name, "x" );
	color_comp = 0;
	color_subcomp = 1;
	iso_number = 5;
	specified_color = 0;
	deform_display_on = 0;
	deform_comp_name[0] = 0;
	deform_comp = 0;
	deform_scale = 0;
	fg_deform_scale = false;
	initial_style = initial_style_none;
	deform_style = deform_style_none;
	initial_line_color[0] = 0;
	initial_line_color[1] = 0;
	initial_line_color[2] = 1;
	deform_line_color[0] = 0;
	deform_line_color[1] = 0;
	deform_line_color[2] = 1;
	strcpy( output_type, "AVS" );

	// for surface_style == surface_style_equivalent ---------

	data_comp_name[0] = 0;
	strcpy( data_subcomp_name, "x" );
	data_comp = 0;
	data_subcomp = 1;
	iso_value = 0;

	// for surface_style == surface_style_user ---------------

	method = method_quadric;
	point[0] = point[1] = point[2] = 0;
	radius = 1;
	length = 1;
	for(int i=0; i<10; i++) coef[i] = 0;

	// rendering parameters(output_type == "BMP") ------------

	CVis_ViewParams::Init();
	isoline_color[0] = isoline_color[1] = isoline_color[2] = 0;
	fg_isoline_color = false;
	boundary_line_on = 0;
}


void CVis_PSR::Write( class CHECData* hecd )
{
	// common ------------------------------------------------

	hecd->WriteLine( "!!common of PSR ----------------------------------");

	WriteVisParam( hecd, "surface_style",  "I", surface_style );
	WriteVisParam( hecd, "display_method", "I", display_method );

	if( color_comp_name[0] != 0 )
		WriteVisParam( hecd, "display_method", "S", color_comp_name );

	if( color_subcomp_name[0] != 0 )
		WriteVisParam( hecd, "color_subcomp_name", "S", color_subcomp_name );

	WriteVisParam( hecd, "color_comp",       "I", color_comp );
	WriteVisParam( hecd, "color_subcomp",    "I", color_subcomp );
	WriteVisParam( hecd, "iso_number",       "I", iso_number );
	WriteVisParam( hecd, "specified_color",  "F", specified_color );
	WriteVisParam( hecd, "deform_display_on","I", deform_display_on );
	if( deform_comp_name[0] != 0 )
		WriteVisParam( hecd, "deform_comp_name", "S", deform_comp_name );
	WriteVisParam( hecd, "deform_comp", "I", deform_comp );
	if( fg_deform_scale )
		WriteVisParam( hecd, "deform_scale", "F", deform_scale );
	if( initial_style != initial_style_none)
		WriteVisParam( hecd, "initial_style", "I", initial_style );
	if( deform_style != deform_style_none)
		WriteVisParam( hecd, "deform_style", "I", deform_style );

	WriteVisPArry( hecd, "initial_line_color", 'F', 3, initial_line_color );
	WriteVisPArry( hecd, "deform_line_color",  'F', 3, deform_line_color );
	WriteVisParam( hecd, "output_type", "S", output_type );

	// for surface_style == surface_style_equivalent ---------

	hecd->WriteLine( "!!for equivalent value surface -------------------");

	if( data_comp_name[0] != 0)
		WriteVisParam( hecd, "data_comp_name", "S", data_comp_name );
	WriteVisParam( hecd, "data_subcomp_name", "S", data_subcomp_name );
	WriteVisParam( hecd, "data_comp",    "I", data_comp );
	WriteVisParam( hecd, "data_subcomp", "I", data_subcomp );
	WriteVisParam( hecd, "iso_value",    "F", iso_value);

	// for surface_style == surface_style_user ---------------

	hecd->WriteLine( "!!for surface of user defined equation -----------");

	WriteVisParam( hecd, "method", "I", method );
	WriteVisPArry( hecd, "point",  'F', 3, point );
	WriteVisParam( hecd, "radius", "F", radius );
	WriteVisParam( hecd, "length", "F", length );
	WriteVisPArry( hecd, "coef",   'F', 10, coef );

	// rendering parameters(output_type == "BMP") ------------
	// parameters defined in CVis_ScreenParams and followings

	hecd->WriteLine( "!!for BMP output ---------------------------------");

	CVis_ViewParams::Write( hecd );
	if( fg_isoline_color )
		WriteVisPArry( hecd, "isoline_color",  'F', 3, isoline_color );
	WriteVisParam( hecd, "boundary_line_on", "I", boundary_line_on );
}



bool CVis_PSR::ReadData( const char* line )
{
	if( CVis_ViewParams::ReadData(line)) return true;

	char buff[256];
	strcpy( buff, line );
	char* header = strtok( buff, " =,\r\n\t" );
	if( header[0] == 0 || header[0] != '!' ) return false;
	tolower( header );

	#define GENERATE_CODE( x, n, fg ) \
		else if(strcmp( #x, header )==0) {\
			int r = ReadVisValue( 0, n, &x );\
			fg = (r>0); \
			return fg; \
		}
	#define GENERATE_CODEA( x, n, fg ) \
		else if(strcmp( #x, header )==0) {\
			int r = ReadVisValue( 0, n, x );\
			fg = (r>0); \
			return fg; \
		}
	bool dummy;
	if(false); // dummy
	GENERATE_CODE( surface_style,       1, dummy )
	GENERATE_CODE( display_method,      1, dummy )
	GENERATE_CODEA( color_comp_name,  100, dummy )
	GENERATE_CODEA( color_subcomp_name, 4, dummy )
	GENERATE_CODE( color_comp,          1, dummy )
	GENERATE_CODE( color_subcomp,       1, dummy )
	GENERATE_CODE( iso_number,          1, dummy )
	GENERATE_CODE( specified_color,     1, dummy )
	GENERATE_CODE( deform_display_on,   1, dummy )
	GENERATE_CODEA( deform_comp_name, 100, dummy )
	GENERATE_CODE( deform_comp,         1, dummy )
	GENERATE_CODE( deform_scale,        1, fg_deform_scale )
	GENERATE_CODE( initial_style,       1, dummy )
	GENERATE_CODE( deform_style,        1, dummy )
	GENERATE_CODEA( initial_line_color, 3, dummy )
	GENERATE_CODEA( deform_line_color,  3, dummy )
	GENERATE_CODEA( output_type,        3, dummy )

	// for surface_style == surface_style_equivalent ---------

	GENERATE_CODEA( data_comp_name,   100, dummy )
	GENERATE_CODEA( data_subcomp_name,  4, dummy )
	GENERATE_CODE( data_comp,           1, dummy )
	GENERATE_CODE( data_subcomp,       1, dummy )
	GENERATE_CODE( iso_value,           1, dummy )

	// for surface_style == surface_style_user ---------------

	GENERATE_CODE( method,              1, dummy )
	GENERATE_CODEA( point,              3, dummy )
	GENERATE_CODE( radius,              1, dummy )
	GENERATE_CODE( length,              1, dummy )
	GENERATE_CODEA( coef,              10, dummy )

	// rendering parameters(output_type == "BMP") ------------
	// parameters defined in CVis_ScreenParams and followings

	GENERATE_CODEA( isoline_color,      3, fg_isoline_color )
	GENERATE_CODE( boundary_line_on,    1, dummy )

	#undef GENERATE_CODE
	#undef GENERATE_CODEA

	return false;
}



//*****************************************************************************
// CVis_PVR
//*****************************************************************************


CVis_PVR::CVis_PVR()
 : CVis_ViewParams(), fea_point(0)
{
	Init();
}

CVis_PVR::~CVis_PVR()
{
	delete fea_point;
}



void CVis_PVR::Init()
{
	// for dividing ------------------------------------------

	maximum_refinement = 100;
	n_voxel_x = 0;
	fg_n_voxel_x = false;
	n_voxel_y = 0;
	fg_n_voxel_y = false;
	n_voxel_z = 0;
	fg_n_voxel_z = false;
	voxel_filename[0] = 0;
	x_specified_level = 100;
	y_specified_level = 100;
	z_specified_level = 100;

	// for opacity

	transfer_function_type = 1;
	opa_value = 0.02;
	num_of_features = 0;
	delete fea_point;
	fea_point = 0;
	name_lookup[0] = 0;

	// parameters defined in CVis_ViewParams and followings

	CVis_ViewParams::Init();
	histogram_on = 0;
	for(int i=0; i<6; i++ ) display_range[i] = 0;
	fg_display_range = false;
}


void CVis_PVR::Write( class CHECData* hecd )
{
	// for dividing ------------------------------------------

	hecd->WriteLine( "!!for dividing -----------------------------------");

	WriteVisParam( hecd, "maximum_refinement", "I", maximum_refinement );
	if( fg_n_voxel_x )
		WriteVisParam( hecd, "n_voxel_x", "I", n_voxel_x );
	if( fg_n_voxel_y )
		WriteVisParam( hecd, "n_voxel_y", "I", n_voxel_y );
	if( fg_n_voxel_z )
		WriteVisParam( hecd, "n_voxel_z", "I", n_voxel_z );
	if( voxel_filename[0] != 0 )
		WriteVisParam( hecd, "voxel_filename", "S", voxel_filename );
	WriteVisParam( hecd, "x_specified_level", "I", x_specified_level );
	WriteVisParam( hecd, "y_specified_level", "I", y_specified_level );
	WriteVisParam( hecd, "z_specified_level", "I", z_specified_level );

	// for opacity -------------------------------------------

	hecd->WriteLine( "!!for opacity ------------------------------------");

	WriteVisParam( hecd, "transfer_function_type", "I", transfer_function_type );
	WriteVisParam( hecd, "opa_value", "F", opa_value );
	WriteVisParam( hecd, "num_of_features", "I", num_of_features );

	if( num_of_features > 0 && fea_point )
		WriteVisPArry( hecd, "fea_point",  'F', 3 * num_of_features, fea_point );

	if( name_lookup[0] != 0 )
		WriteVisParam( hecd, "name_lookup", "S", name_lookup );

	// parameters defined in CVis_ViewParams and followings

	hecd->WriteLine( "!!view params ------------------------------------");

	CVis_ViewParams::Write( hecd );

	WriteVisParam( hecd, "histogram_on", "I", histogram_on );
	if( fg_display_range )
		WriteVisPArry( hecd, "display_range",  'F', 6, display_range );
}


bool CVis_PVR::ReadData( const char* line )
{
	if( CVis_ViewParams::ReadData(line)) return true;

	char buff[256];
	strcpy( buff, line );
	char* header = strtok( buff, " =,\r\n\t" );
	if( header[0] == 0 || header[0] != '!' ) return false;
	tolower( header );

	#define GENERATE_CODE( x, n, fg ) \
		else if(strcmp( #x, header )==0) {\
			int r = ReadVisValue( 0, n, &x );\
			fg = (r>0); \
			return fg; \
		}
	#define GENERATE_CODEA( x, n, fg ) \
		else if(strcmp( #x, header )==0) {\
			int r = ReadVisValue( 0, n, x );\
			fg = (r>0); \
			return fg; \
		}

	bool dummy;
	if(false); // dummy

	// for dividing ------------------------------------------

	GENERATE_CODE( maximum_refinement,  1, dummy )
	GENERATE_CODE( n_voxel_x,  1, fg_n_voxel_x )
	GENERATE_CODE( n_voxel_y,  1, fg_n_voxel_y )
	GENERATE_CODE( n_voxel_z,  1, fg_n_voxel_z )
	GENERATE_CODEA( voxel_filename, 100, dummy )
	GENERATE_CODE( x_specified_level, 1, dummy )
	GENERATE_CODE( y_specified_level, 1, dummy )
	GENERATE_CODE( z_specified_level, 1, dummy )

	// for opacity

	GENERATE_CODE( transfer_function_type, 1, dummy )
	GENERATE_CODE( opa_value,              1, dummy )
	GENERATE_CODE( num_of_features,        1, dummy )
	GENERATE_CODEA( name_lookup,         100, dummy )

	// parameters defined in CVis_ViewParams and followings

	GENERATE_CODE( histogram_on,           1, dummy )
	GENERATE_CODEA( display_range,         6, fg_display_range )

	#undef GENERATE_CODE
	#undef GENERATE_CODEA

	#define GENERATE_CODE( x, n, type, fmt ) \
		else if(strcmp( #x, header )==0) {\
			delete[] x; x = 0; \
			x = ReadVisPArry<type>( 0, n, fmt );\
			return x!=0; \
		}

	if( num_of_features <= 0 ) return false;

	int f_data_size = 3 * num_of_features;
	if(false);
	GENERATE_CODE( fea_point, f_data_size, double, "%lf" );
	#undef GENERATE_CODE

	return false;
}



//-----------------------------------------------------------------------------
// CHECDB_Visual
//-----------------------------------------------------------------------------


CHECDB_Visual::CHECDB_Visual()
 : CHECDataBlock(HECDB_VISUAL), psr(), pvr(0)
{
	visual_start_step = -1;
	visual_interval_step = 1;
	visual_end_step = -1;
 	surface_num = 0;
}


CHECDB_Visual::~CHECDB_Visual()
{
	Clear();
}


void CHECDB_Visual::Clear()
{
	visual_start_step = -1;
	visual_interval_step = 1;
	visual_end_step = -1;

 	surface_num = 0;

	vector<CVis_PSR*>::iterator is;
	for(is = psr.begin(); is != psr.end(); is++) {
		delete *is;
	}
	psr.clear();

	delete pvr;
	pvr = 0;
}


void CHECDB_Visual::Write( class CHECData* hecd )
{
	surface_num = psr.size();

	if( surface_num > 0 ){
		hecd->WriteHeader( "!VISUAL", "SIII",
			"method",              "PSR",
			"visual_start_step",    visual_start_step,
			"visual_interval_step", visual_interval_step,
			"visual_end_step",      visual_end_step
		);
		CVis_ViewParams::WriteVisParam( hecd, "surface_num", "I", surface_num );
		vector<CVis_PSR*>::iterator is;
		for(is = psr.begin(); is != psr.end(); is++ ) {
			hecd->WriteLine( "!SURFACE" );
			(*is)->Write( hecd );
		}
	}

	if( pvr ){
		hecd->WriteHeader( "!VISUAL", "SIII",
			"method",               "PVR",
			"visual_start_step",    visual_start_step,
			"visual_interval_step", visual_interval_step,
			"visual_end_step",      visual_end_step
		);
		pvr->Write( hecd );
	}
}


bool CHECDB_Visual::Read( CHECData* hecd, char* header_line )
{
	char s[256];
	char method_s[256];
	int rcode[10];

	if(!hecd->ParseHeader( header_line, rcode, "SIII",
		"method",              s,
		"visual_start_step",    &visual_start_step,
		"visual_interval_step", &visual_interval_step,
		"visual_end_step",      &visual_end_step
	)) return false;

	cleanup_token( s, method_s );
	toupper( method_s );
	bool fg_psr = ( strcmp( "PSR", method_s ) == 0 );

	if( fg_psr ) {
		if(! hecd->ReadParameter( rcode, "I", "surface_num", &surface_num )) return false;
		for(int i=0; i<surface_num; i++ ){
			if(! hecd->ReadLine(s)) return false;
			cleanup_token(s);
			toupper(s);
			if( strcmp( "!SURFACE", s ) != 0 ) return false;
			CVis_PSR* vis_psr = new CVis_PSR();
			if(!vis_psr->Read( hecd )) {
				delete vis_psr;
				return false;
			}
			psr.push_back( vis_psr );
		}
	} else {
		delete pvr;
		pvr = new CVis_PVR();
		if(!pvr->Read( hecd )) {
			delete pvr;
			pvr = 0;
			return false;
		}
	}

	return true;
}


