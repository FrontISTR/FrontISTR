/*=====================================================================*
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : I/O and Utility                                   !
!                                                                      !
!            Written by Noboru Imai (Univ. of Tokyo)                   !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
 *=====================================================================*/

#define fstr_ctrl_util_MAIN

#include "fstr_ctrl_util.h"
#include "hecmw_malloc.h"


#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif


#define buffsize       256
#define STR_SIZE       buffsize
#define FILE_NAME_SIZE 512

char err_msg[buffsize];

void c_fstr_ctrl_get_err_msg(char* buff)
{
	strcpy( buff, err_msg );
}


static void set_param_err_msg( fstr_ctrl_data* ctrl, const char* param_name, const char* msg );
static void set_data_err_msg( fstr_ctrl_data* ctrl, int line, int data_no, const char* msg );
static void set_record_data_line_err_msg( fstr_ctrl_data* ctrl, int r );
static char* gettoken( const char* line );
static void strcpy_f2c( char* dest, const char* src, int len);
static void strcpy_c2f( char* dest, int len, const char* src);
static char* remove_header_space( char* token );
static int Strncmpi( const char* s1, const char* s2, int len);
static int Strcmpi( const char* s1, const char* s2 );
/* static void Strupr( char* s ); */
static void remove_cr( char* s );
static void format_conv( const char* format, char* fmt, int* array_size );

/* JP-0 */
/* JP-1 */

int fg_fortran_get_data_array_v = 0;

/* ================================================================================= */
/* STATIC FUNCTIONS                                                                  */
/* ================================================================================= */


static void fstr_ctrl_data_init( fstr_ctrl_data* ctrl )
{
	ctrl->rec = NULL;
	ctrl->rec_n = 0;
	ctrl->header_pos = NULL;
	ctrl->header_n = 0;
	ctrl->data_line_n = NULL;
	ctrl->current_header_index = 0;

	err_msg[0] = 0;
}


static void fstr_ctrl_data_finalize( fstr_ctrl_data* ctrl )
{
	int i;

	for( i=0; i<ctrl->rec_n;i++ ){
		HECMW_free(ctrl->rec[i].line);
	}
	HECMW_free(ctrl->rec);
	HECMW_free(ctrl->header_pos);
	HECMW_free(ctrl->data_line_n);
}


int fstr_ctrl_tmp_rewind( fstr_ctrl_data* ctrl )
{
	ctrl->current_header_index = 0;
	return 0;
}

int fstr_ctrl_rewind__( int* ctrl )
{
	return fstr_ctrl_tmp_rewind( ctrl_list[*ctrl] );
}

int fstr_ctrl_rewind( int* ctrl ) { return fstr_ctrl_rewind__( ctrl ); }
int fstr_ctrl_rewind_( int* ctrl ) { return fstr_ctrl_rewind__( ctrl ); }
int FSTR_CTRL_REWIND( int* ctrl ) { return fstr_ctrl_rewind__( ctrl ); }
int FSTR_CTRL_REWIND_( int* ctrl ) { return fstr_ctrl_rewind__( ctrl ); }
int FSTR_CTRL_REWIND__( int* ctrl ) { return fstr_ctrl_rewind__( ctrl ); }


/*-----------------------------------------------------------------------------------*/


int line_no;

static int is_comment_or_blank_line( char* buff )
{
	char* p = buff;

	if(  buff[0] == '\n' || buff[0] == 0 || buff[0] == '#' || (buff[0]=='!' && buff[1] == '!'))
		return TRUE;

	while( *p ) {
		if( *p != ' ' && *p !='\n' )
			return FALSE;
		p++;
	}
	return TRUE;
}

/*-----------------------------------------------------------------------------------*/


static int count_line_and_header_number( const char* fname, int* line_n, int* header_n )
{
	FILE* fp;
	char buff[buffsize];
	int line_no, L, N;

	fp = fopen( fname, "r");
	if(!fp){
		return -1;
	}

	line_no = 0;
	L = N = 0;

	while( !feof(fp)){
		line_no++;
		if( fgets( buff, buffsize-1, fp) == NULL ) break;
		if( is_comment_or_blank_line( buff ))
			continue;
		if( buff[0] == '!' ) {
			N++;
		}
		L++;
	}
	fclose(fp);

	*line_n = L;
	*header_n = N;
	return 0;
}


/*-----------------------------------------------------------------------------------*/



static int set_fstr_ctrl_data( const char* fname, fstr_ctrl_data* ctrl )
{
	FILE* fp;
	char buff[buffsize];
	int line_no;
	int header_count;
	int rec_count;

	fp = fopen( fname, "r");

	if(!fp){
		return -1;
	}

	line_no = 0;
	rec_count = 0;
	header_count = -1;

	while( fgets( buff, buffsize-1, fp) ){
		line_no++;

		if( is_comment_or_blank_line( buff )){
			continue;
		}

		if( buff[0] == '!' ) {
			header_count++;
			ctrl->header_pos[ header_count] = rec_count;
			ctrl->data_line_n[ header_count ] = 0;
		} else {
			if( header_count>=0 ){
				ctrl->data_line_n[ header_count ] ++;
			}
		}
		ctrl->rec[rec_count].line_no = line_no;
		ctrl->rec[rec_count].line = HECMW_malloc( sizeof(char) * (strlen( buff )+1));

		if( ctrl->rec[rec_count].line == NULL ) {
			printf("Not enough memory\n");
			exit(-1);
		}
		/* Strupr( buff ); */
		remove_cr( buff );
		strcpy( ctrl->rec[rec_count].line, buff );

		rec_count++;
	}

	fclose(fp);
	return 0;
}

/*-----------------------------------------------------------------------------------*/

static int create_fstr_ctrl_data( const char* fname, fstr_ctrl_data* ctrl )
{
	int i;
	int line_n, header_n;

	fstr_ctrl_data_init( ctrl );
	if( count_line_and_header_number( fname, &line_n, &header_n )) return -1;

	ctrl->rec_n = line_n;
	ctrl->rec = (ctrl_rec*)HECMW_malloc( sizeof(ctrl_rec) * line_n );
	for(i=0; i<line_n; i++) {
		ctrl->rec[i].line = NULL;
		ctrl->rec[i].line_no = 0;
	}

	ctrl->header_n = header_n;
	ctrl->header_pos = (int*)HECMW_malloc( sizeof(int) * header_n );
	ctrl->data_line_n = (int*)HECMW_malloc( sizeof(int) * header_n );
	for(i=0; i<header_n; i++) {
		ctrl->header_pos[i] = 0;
		ctrl->data_line_n[i] = 0;
	}

	if( set_fstr_ctrl_data( fname, ctrl )) {
		fstr_ctrl_data_finalize( ctrl );
		return -1;
	}

	return 0;
}


/*===================================================================================*/
/* PUBLIC FUNCTIONS                                                                  */
/*===================================================================================*/


/* return fstr_ctrl_data*  --- NULL: error */


fstr_ctrl_data* c_fstr_ctrl_open( const char* filename )
{
	fstr_ctrl_data* ctrl;

	ctrl = (fstr_ctrl_data*)HECMW_malloc( sizeof( fstr_ctrl_data ) );
	if(!ctrl) return NULL;

#if 0  /* all PE reading fstr contorl */
/* #ifndef HECMW_SERIAL */
	{
		int i;
		int myrank = HECMW_comm_get_rank();
		HECMW_Comm comm = HECMW_comm_get_comm();
		/* file ----------------------------------------------------------- */
		int err;
		if( myrank == 0 ) {
			err = create_fstr_ctrl_data( filename, ctrl );
		}
		HECMW_Bcast( &err, 1, HECMW_INT, 0, comm );
		if( err ) {
			HECMW_free( ctrl );
			return NULL;
		}

		/* rec ------------------------------------------------------------ */

		HECMW_Bcast( &ctrl->rec_n,    1, HECMW_INT, 0, comm );
		if( myrank != 0 ){
			ctrl->rec = HECMW_malloc( sizeof( ctrl_rec ) * ctrl->rec_n );
		}
		for( i=0; i<ctrl->rec_n; i++ ){
			int line_size;
			if( myrank == 0 )
				line_size = strlen( ctrl->rec[i].line ) + 1;
			HECMW_Bcast( &line_size,  1, HECMW_INT, 0, comm );
			if( myrank != 0 ) {
				ctrl->rec[i].line = HECMW_malloc( sizeof(char) * line_size );
			}
			HECMW_Bcast( ctrl->rec[i].line,    line_size, HECMW_CHAR, 0, comm );
			HECMW_Bcast( &ctrl->rec[i].line_no,        1, HECMW_INT,  0, comm );
		}

		/* header --------------------------------------------------------- */

		HECMW_Bcast( &ctrl->header_n, 1, HECMW_INT, 0, comm );
		if( myrank != 0 ){
			ctrl->header_pos = HECMW_malloc( sizeof( int ) * ctrl->header_n );
		}
		HECMW_Bcast( ctrl->header_pos, ctrl->header_n, HECMW_INT, 0, comm );

		/* data line  ----------------------------------------------------- */

		if( myrank != 0 ){
			ctrl->data_line_n = HECMW_malloc( sizeof( int ) * ctrl->header_n );
		}
		HECMW_Bcast( ctrl->data_line_n, ctrl->header_n, HECMW_INT, 0, comm );
	}
#else
	if( create_fstr_ctrl_data( filename, ctrl )) {
		HECMW_free( ctrl );
		return NULL;
	}
#endif
	return ctrl;
}


/*-----------------------------------------------------------------------------------*/
/* JP-2 */

int c_fstr_ctrl_get_rec_number( fstr_ctrl_data* ctrl )
{
	if( !ctrl ) return -1;

	return ctrl->rec_n;
}

/*-----------------------------------------------------------------------------------*/
/* JP-3 */

int c_fstr_ctrl_get_line( fstr_ctrl_data* ctrl, int rec_no, char* buff )
{
	if( !ctrl ) return -1;
	if( rec_no < 0 || rec_no >= ctrl->rec_n ) return -1;

	strcpy( buff, ctrl->rec[ rec_no ].line );
	return 0;
}

/*-----------------------------------------------------------------------------------*/
/* JP-4 */

int c_fstr_ctrl_seek_header( fstr_ctrl_data* ctrl, const char* header_name )
{
	int i;
	static char* h_name = NULL;
	int header_name_len;
	int start_index;

	if( !ctrl ) return FALSE;

	if( header_name == NULL || header_name[0]==0 ){
		start_index = ctrl->current_header_index + 1;
	} else {
		start_index = 0;
		h_name = (char*)header_name;
	}

	if( h_name == NULL )
		return FALSE;

	header_name_len = strlen( h_name );

	for( i = start_index; i< ctrl->header_n; i++){
		int hp = ctrl->header_pos[i];
		char* header = ctrl->rec[hp].line;
		if( Strncmpi( header, h_name, header_name_len) == 0) {
			ctrl->current_header_index = i;
			return TRUE;
		}
	}
	return FALSE;
}

/*-----------------------------------------------------------------------------------*/
/* JP-5 */

int c_fstr_ctrl_seek_next_header( fstr_ctrl_data* ctrl )
{

	if( !ctrl ) return FALSE;

	ctrl->current_header_index++;
	if( ctrl->header_n <= ctrl->current_header_index )
		return FALSE;
	else
		return TRUE;
}

/*-----------------------------------------------------------------------------------*/
/* JP-6 */

int c_fstr_ctrl_get_current_header_name( fstr_ctrl_data* ctrl, char* header_name )
{
	int hp;
	char* line_p;
	char* header_p;


	if( !ctrl ) return -1;
	hp = ctrl->header_pos[ ctrl->current_header_index ];

	line_p = ctrl->rec[hp].line;
	header_p = header_name;

	while( *line_p && ( *line_p != ',' ) && ( *line_p != ' ') && ( *line_p != '\n') && ( *line_p != '\r')){
		*header_p = (char)toupper(*line_p);
		line_p++;
		header_p++;
	}
	*header_p = 0;

	return 0;
}

/*-----------------------------------------------------------------------------------*/
/* JP-7 */

int c_fstr_ctrl_get_current_header_line_no( fstr_ctrl_data* ctrl )
{
	int index, hp;

	if( !ctrl ) return -1;

	index = ctrl->current_header_index;
	hp = ctrl->header_pos[index];
	return ctrl->rec[hp].line_no;
}

/*-----------------------------------------------------------------------------------*/
/* JP-8 */

int c_fstr_ctrl_get_current_header_pos( fstr_ctrl_data* ctrl )
{
	int index, hp;

	if( !ctrl ) return -1;

	index = ctrl->current_header_index;
	hp = ctrl->header_pos[index];
	return hp;
}




/*-----------------------------------------------------------------------------------*/
/* JP-9 */


static int param_value_check( const char* value_list, char* val, int *index )
{
	char vlist[buffsize];
	char* token;
	int fg, i, n;

	if( value_list == NULL || value_list[0] == 0 || value_list[0] == '#'){
		if(index)
			*index = 0;
		return 0;
	}

	/* JP-10 */
	fg = 0;
	n = strlen( value_list );
	for( i=0; i<n; i++) {
		if(value_list[i] != ' ') {
			fg = 1;
			break;
		}
	}
	if( !fg ){
		if(index)
			*index = 0;
		return 0;
	}

	strcpy( vlist, value_list );

	i = 1;
	token = strtok( vlist, " ,");
	while(token){
		if( Strcmpi( token, val ) == 0) {
			if(index)
				*index = i;
			return 0;
		}
		token = strtok( NULL, ", ");
		i++;
	}
	return -1;
}

static int param_value_convert( int type, char* token, void* val )
{
	char fmt[5];
	int fmt_i;
	int r;

	fmt_i = 0;
	fmt[fmt_i] = '%'; fmt_i++;
	type = toupper(type);
	switch( (char)type ) {
	case 'I':
		fmt[fmt_i] = 'd'; fmt_i++;
		break;
	case 'C':
	case 'S':
	case (int)'P':
		fmt[fmt_i] = 's'; fmt_i++;
		break;
	case 'R':
		fmt[fmt_i] = 'l'; fmt_i++;
		fmt[fmt_i] = 'f'; fmt_i++;
		break;
	default:
		return FSTR_CTRL_RCODE_PARAM_TYPE_ERROR;
	}
	fmt[fmt_i] = 0;
	r = sscanf( token, fmt, val );
	if( r != 1 )
		return FSTR_CTRL_RCODE_PARAM_TYPE_ERROR;

	return FSTR_CTRL_RCODE_PARAM_SUCCESS;
}



int c_fstr_ctrl_get_param( fstr_ctrl_data* ctrl, const char* param_name, const char* value_list, char type, void* val )
{
	int h_index;
	int h_pos;
	char header[buffsize];

	type = type & 0xff; /* bug fix for compiler=pgi */
	if( !ctrl ) return -1;

	h_index = ctrl->current_header_index;
	if( h_index < 0 ) return -1;

	h_pos = ctrl->header_pos[h_index];
	strcpy( header, ctrl->rec[h_pos].line);

	{
		char* token;
		char* param_pos;
		char* eq_pos;
		char* val_pos;
		char* p;
		int param_name_len;
		int r;
		int index;

		param_name_len = strlen(param_name);

		strtok( header, ",\n"); /* remove !header name  */
		token = strtok( NULL, ",\n");
		while( token ) {
			param_pos = remove_header_space( token );
			if( Strncmpi(param_pos, param_name, param_name_len)==0 ){
				if( type == 'E' || type == 'e' ){
					*((int*)val) = 1;
					return FSTR_CTRL_RCODE_PARAM_VALUE_NOTHING;
				}
				eq_pos = strstr( param_pos, "=");
				if( eq_pos ){
					val_pos = eq_pos + 1;
					val_pos = remove_header_space( val_pos );
					p = val_pos;
					while( *p ) {
						if( *p == ','){
							*p = 0;
							break;
						}
						p++;
					}

					if( param_value_check( value_list, val_pos, &index)){
						return FSTR_CTRL_RCODE_PARAM_RANGE_ERROR;
					}
					if( type == 'P' || type == 'p' ){
						*((int*)val) = index;
						return FSTR_CTRL_RCODE_PARAM_SUCCESS;
					} else {
						r = param_value_convert( type, val_pos, val);
						return r;
					}
				} else {
					return FSTR_CTRL_RCODE_PARAM_VALUE_NOTHING;
				}
			}
			token = strtok( NULL, ",\n");
		}
	}

	if( type == 'E' || type == 'e' ){
		*((int*)val) = 0;
	}
	return FSTR_CTRL_RCODE_PARAM_NOTHING;
}

/*-----------------------------------------------------------------------------------*/
/* JP-11 */


static int rcode_of_get_param = 0;

int c_fstr_ctrl_get_param_ex( fstr_ctrl_data* ctrl,
		const char* param_name, const char* value_list, int necessity, char type, void* val )
{
	char s[buffsize];
	int rcode = c_fstr_ctrl_get_param(ctrl, param_name, value_list, type, val );
	rcode_of_get_param  = rcode;

	switch( rcode ){
	case FSTR_CTRL_RCODE_PARAM_SUCCESS:
		return FSTR_CTRL_RCODE_PARAM_SUCCESS;
	case FSTR_CTRL_RCODE_PARAM_NOTHING: /* nothing parameter */
		if( necessity ) {
			set_param_err_msg( ctrl, param_name, " is required");
		} else
			return FSTR_CTRL_RCODE_PARAM_SUCCESS;
		break;
	case FSTR_CTRL_RCODE_PARAM_VALUE_NOTHING: /* nothing value of parameter */
		if( necessity ) {
			set_param_err_msg( ctrl, param_name, ": value is required");
		} else
			return FSTR_CTRL_RCODE_PARAM_SUCCESS;
		break;
	case FSTR_CTRL_RCODE_PARAM_TYPE_ERROR: /* type change error */
		set_param_err_msg( ctrl, param_name, ": type conversion fail");
		break;
	case FSTR_CTRL_RCODE_PARAM_RANGE_ERROR: /* range error */
		sprintf( s, ": type range fail(%s)", value_list );
		set_param_err_msg( ctrl, param_name, s );
		break;
	default:
		assert(0);
	}

	return rcode;
}



/*-----------------------------------------------------------------------------------*/
/* JP-12 */

int c_fstr_ctrl_get_data_line_n( fstr_ctrl_data* ctrl )
{
	int h_index;

	if( !ctrl ) return -1;

	h_index = ctrl->current_header_index;
	if( h_index < 0 ) return -1;

	return  ctrl->data_line_n[h_index];
}

/*-----------------------------------------------------------------------------------*/

/* JP-13 */

int c_fstr_ctrl_copy_data_line( fstr_ctrl_data* ctrl, int line_no, char* data_line )
{
	int data_line_n;
	int h_index;
	int data_pos;

	data_line_n = c_fstr_ctrl_get_data_line_n( ctrl );
	if( data_line_n <= 0 || data_line_n < line_no )
		return -1;

	h_index = ctrl->current_header_index;
	data_pos = ctrl->header_pos[ h_index ] + line_no;

	strcpy( data_line, ctrl->rec[data_pos].line);
	return 0;
}


/*-----------------------------------------------------------------------------------*/

/* JP-14 */

int c_fstr_ctrl_get_data_n_in_line( fstr_ctrl_data* ctrl, int line_no, const char* delim )
{
	char data_line[ buffsize ];
	char* token;
	int err;
	int counter;

	err = c_fstr_ctrl_copy_data_line( ctrl, line_no, data_line );
	if( err ) return -1;

	counter = 0;
	token = strtok( data_line, delim );
	while( token ){
		counter++;
		token = strtok( NULL, delim );
	}

	return counter;
}


/*-----------------------------------------------------------------------------------*/


/* JP-15 */

int error_pos = -1;

int c_fstr_ctrl_get_data_error_pos(void)
{
	return error_pos;
}


int error_line = -1;

int c_fstr_ctrl_get_data_error_line(void)
{
	return error_line;
}


/* JP-16 */


int c_fstr_ctrl_get_data( fstr_ctrl_data* ctrl, int line_no, const char* format, ... )
{
	va_list va;
	int r;

	va_start(va, format);
	r = c_fstr_ctrl_get_data_v( ctrl, line_no, format, va );
	va_end(va);

	return r;
}



int c_fstr_ctrl_get_data_ex( fstr_ctrl_data* ctrl, int line_no, const char* format, ... )
{
	va_list va;
	int r;

	va_start(va, format);
	r = c_fstr_ctrl_get_data_v( ctrl, line_no, format, va );
	va_end(va);

	if( r != 0 ) {
		set_record_data_line_err_msg( ctrl, r );
		return -1;
	}

	return 0;
}



int c_fstr_ctrl_get_data_v( fstr_ctrl_data* ctrl, int line_no, const char* format, va_list va )
{
	char data_line[ buffsize ];
	int necessary;
	int type;
	void* val_p;
	char* token;
	int err;
	int counter;
	int len;
	char* fmt;
	char fmt_integer[] = "%d";
	char fmt_double[]  = "%lf";
	char fmt_string[]  = "%s";
	char fmt_char[]    = "%c";
        char format_c[buffsize];
        int array_size[buffsize];
        char buff[buffsize*2];

        format_conv(format, format_c, array_size);

	error_pos = -1;

	err = c_fstr_ctrl_copy_data_line( ctrl, line_no, data_line );
	if( err ) {
		int i=0;
		while( format_c[i] != 0 ){
			if( isupper(format_c[i]) ){
				return FSTR_CTRL_RCODE_DATA_NOTHING;
			}
			i++;
		}
		return 0;
	}

	len = strlen( format_c );
	counter = 0;

	if( data_line[0] == '!' ){
		for( ; counter<len; counter++ ) {
			if( isupper(format_c[counter])){
				return counter+1; 
			}
		}
		return 0;
	}


	token = gettoken( data_line );
	while( token && counter<len){
		error_pos = counter+1;
		necessary = isupper(format_c[counter]);
		type = toupper(format_c[counter]) & 0xff;
		switch( (char)type ){
		case 'I':
			fmt = fmt_integer;
			break;
		case 'R':
			fmt = fmt_double;
			break;
		case 'S':
			fmt = fmt_string;
			break;
		case 'C':
			fmt = fmt_char;
			break;
		default:
			return FSTR_CTRL_RCODE_DATA_ERROR;
		}
		val_p = va_arg( va, void* );
		if( token[0] != 0 ) {
			err = sscanf( token, fmt, val_p );
			if( err != 1 && necessary ) {
				return FSTR_CTRL_RCODE_DATA_TYPE_ERROR;
			}
                        if( (char)type == 'S' ) {
                                strcpy( buff, (char*)val_p );
                                strcpy_c2f( (char*)val_p, array_size[counter], buff );
                        }
		} else if( necessary ) {
			return FSTR_CTRL_RCODE_DATA_NOTHING;
		}
		token = gettoken( NULL );
		counter++;
	}

	for( ; counter<len; counter++ ) {
		if( isupper(format[counter])){
			return FSTR_CTRL_RCODE_DATA_NOTHING;
		}
	}

	return FSTR_CTRL_RCODE_DATA_SUCCESS;
}





int c_fstr_ctrl_get_data_array_v( fstr_ctrl_data* ctrl, const char* format, va_list va )
{
	#define MAX_DATA_ARRAY_NUMBER 20

	int line_n;
	int i,j;
	char fmt[buffsize];
	int array_size[buffsize];
	char *param[MAX_DATA_ARRAY_NUMBER];
	size_t param_size[MAX_DATA_ARRAY_NUMBER];
	char* p;
	int column_n;
	int type;
	int r;
	char buff[buffsize*2];

	line_n = c_fstr_ctrl_get_data_line_n( ctrl );
	if( line_n < 0 )
		return FSTR_CTRL_RCODE_DATA_LINE_NOTHING;

	format_conv( format, fmt, array_size );
	p = fmt;
	i = 0;
	while(*p){
		if( i>= MAX_DATA_ARRAY_NUMBER ){
			assert(0);
		}
		param[i] = va_arg( va, void* );
		type = toupper( *p ) & 0xff;
		switch((char)type){
		case 'I':
			param_size[i] = sizeof(int);
			break;
		case 'R':
			param_size[i] = sizeof(double);
			break;
		case 'S':
			param_size[i] = sizeof(char) * array_size[i];
			break;
		case 'C':
			param_size[i] = sizeof(char);
			break;
		default:
			assert(0);
		}
		
		p++;
		i++;
	}
	column_n = i;
	for(i=1; i<=line_n; i++){
		r = c_fstr_ctrl_get_data( ctrl, i, fmt,
			param[0],  param[1],  param[2],  param[3],  param[4],
			param[5],  param[6],  param[7],  param[8],  param[9],
			param[10], param[11], param[12], param[13], param[14],
			param[15], param[16], param[17], param[18], param[19] );
		if( r != 0 ) {
			error_line = i;
			return r;
		}
		for(j=0; j<column_n; j++){
			if( fg_fortran_get_data_array_v && toupper(fmt[j]) == 'S' ){
				strcpy( buff, (char*)param[j] );
				strcpy_c2f( (char*)param[j], array_size[j], buff );
			}
			param[j] += param_size[j];
		}
	}

	return FSTR_CTRL_RCODE_DATA_SUCCESS;

	#undef MAX_DATA_ARRAY_NUMBER
}


int c_fstr_ctrl_get_data_array( fstr_ctrl_data* ctrl, const char* format, ... )
{
	va_list va;
	int rcode;

	va_start( va, format );
	rcode = c_fstr_ctrl_get_data_array_v( ctrl, format, va );
	va_end(va);

	return rcode;
}



int c_fstr_ctrl_get_data_array_ex( fstr_ctrl_data* ctrl, const char* format, ... )
{
	va_list va;
	int r;

	va_start( va, format );
	r = c_fstr_ctrl_get_data_array_v( ctrl, format, va );
	va_end(va);

	if( r != 0 ) {
		set_record_data_line_err_msg( ctrl, r);
		return -1;
	}

	return 0;
}


/*-----------------------------------------------------------------------------------*/

/* JP-17 */


int c_fstr_ctrl_close( fstr_ctrl_data* ctrl )
{
	fstr_ctrl_data_finalize( ctrl );
	return 0;
}


/*-----------------------------------------------------------------------------------*/


void c_fstr_ctrl_dump(  fstr_ctrl_data* ctrl )
{
	int i;

	printf("header pos: ");
	for(i=0; i<ctrl->header_n; i++){
		printf("%d,", ctrl->header_pos[i]);
	}
	printf("\n");

	for(i=0; i<ctrl->rec_n; i++){
		printf("%3d: %s", ctrl->rec[i].line_no, ctrl->rec[i].line);
	}
}



/* ==================================================================================*/
/* FORTRAN INTERFACE                                                                 */
/* ==================================================================================*/

void fstr_ctrl_get_err_msg( char* f_buff, int* len )
{
	strcpy_c2f( f_buff, *len, err_msg);
}

void fstr_ctrl_get_err_msg_( char* f_buff, int* len )
{
	strcpy_c2f( f_buff, *len, err_msg);
}

void fstr_ctrl_get_err_msg__( char* f_buff, int* len )
{
	strcpy_c2f( f_buff, *len, err_msg);
}

void FSTR_CTRL_GET_ERR_MSG( char* f_buff, int* len )
{
	strcpy_c2f( f_buff, *len, err_msg);
}

void FSTR_CTRL_GET_ERR_MSG_( char* f_buff, int* len )
{
	strcpy_c2f( f_buff, *len, err_msg);
}

void FSTR_CTRL_GET_ERR_MSG__( char* f_buff, int* len )
{
	strcpy_c2f( f_buff, *len, err_msg);
}

/*-----------------------------------------------------------------------------------*/


/* JP-18 */

int is_first = 1; /* JP-19 */


int fstr_ctrl_open( char* filename )
{
	int i;
	int index;
	char fname[FILE_NAME_SIZE];
	strcpy_f2c( fname, filename, FILE_NAME_SIZE );

	if( is_first ) {
		for(i =0; i<ctrl_list_size; i++) {
			ctrl_list[i] = NULL;
		}
		index = 0;
		is_first = 0;
	} else {
		index = -1;
		for(i =0; i<ctrl_list_size; i++) {
			if( ctrl_list[i] == NULL ) {
				index = i;
				break;
			}
		}
		if( index < 0 ) return -1;
	}
	ctrl_list[index] = c_fstr_ctrl_open( fname );
	if( ctrl_list[index] == NULL )
		return -1;
	else
		return index;
}

int fstr_ctrl_open_( char* filename ) { return fstr_ctrl_open( filename ); }
int fstr_ctrl_open__( char* filename ) { return fstr_ctrl_open( filename ); }
int FSTR_CTRL_OPEN( char* filename ) { return fstr_ctrl_open( filename ); }
int FSTR_CTRL_OPEN_( char* filename ) { return fstr_ctrl_open( filename ); }
int FSTR_CTRL_OPEN__( char* filename ) { return fstr_ctrl_open( filename ); }


/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_rec_number( int* ctrl ){
	return c_fstr_ctrl_get_rec_number( ctrl_list[*ctrl] );
}

int fstr_ctrl_get_rec_number_( int* ctrl ){
	return c_fstr_ctrl_get_rec_number( ctrl_list[*ctrl] );
}

int fstr_ctrl_get_rec_number__( int* ctrl ){
	return c_fstr_ctrl_get_rec_number( ctrl_list[*ctrl] );
}

int FSTR_CTRL_GET_REC_NUMBER( int* ctrl ){
	return c_fstr_ctrl_get_rec_number( ctrl_list[*ctrl] );
}

int FSTR_CTRL_GET_REC_NUMBER_( int* ctrl ){
	return c_fstr_ctrl_get_rec_number( ctrl_list[*ctrl] );
}

int FSTR_CTRL_GET_REC_NUMBER__( int* ctrl ){
	return c_fstr_ctrl_get_rec_number( ctrl_list[*ctrl] );
}

/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_line( int* ctrl, int* rec_no, char* buff, int* buff_size )
{
	char c_buff[STR_SIZE];

	if( c_fstr_ctrl_get_line( ctrl_list[*ctrl], *rec_no, c_buff )) return -1;

	strcpy_c2f( buff, *buff_size, c_buff);
	return 0;
}


int fstr_ctrl_get_line_( int* ctrl, int* rec_no, char* buff, int* buff_size ) {
	return fstr_ctrl_get_line( ctrl, rec_no, buff, buff_size );
}

int fstr_ctrl_get_line__( int* ctrl, int* rec_no, char* buff, int* buff_size ) {
	return fstr_ctrl_get_line( ctrl, rec_no, buff, buff_size );
}

int FSTR_CTRL_GET_LINE( int* ctrl, int* rec_no, char* buff, int* buff_size ) {
	return fstr_ctrl_get_line( ctrl, rec_no, buff, buff_size );
}

int FSTR_CTRL_GET_LINE_( int* ctrl, int* rec_no, char* buff, int* buff_size ) {
	return fstr_ctrl_get_line( ctrl, rec_no, buff, buff_size );
}

int FSTR_CTRL_GET_LINE__( int* ctrl, int* rec_no, char* buff, int* buff_size ) {
	return fstr_ctrl_get_line( ctrl, rec_no, buff, buff_size );
}

/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_seek_header( int* ctrl, const char* header_name )
{
	char name[STR_SIZE];
	strcpy_f2c( name, header_name, STR_SIZE );

	if( name[0] == 0 )
		return c_fstr_ctrl_seek_header( ctrl_list[*ctrl], NULL );
	else
		return c_fstr_ctrl_seek_header( ctrl_list[*ctrl], name );
}

int fstr_ctrl_seek_header_( int* ctrl, const char* header_name ) {
	return fstr_ctrl_seek_header( ctrl, header_name );
}

int fstr_ctrl_seek_header__( int* ctrl, const char* header_name ) {
	return fstr_ctrl_seek_header( ctrl, header_name );
}

int FSTR_CTRL_SEEK_HEADER( int* ctrl, const char* header_name ) {
	return fstr_ctrl_seek_header( ctrl, header_name );
}

int FSTR_CTRL_SEEK_HEADER_( int* ctrl, const char* header_name ) {
	return fstr_ctrl_seek_header( ctrl, header_name );
}

int FSTR_CTRL_SEEK_HEADER__( int* ctrl, const char* header_name ) {
	return fstr_ctrl_seek_header( ctrl, header_name );
}


/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_seek_next_header( int* ctrl ) {
	return c_fstr_ctrl_seek_next_header( ctrl_list[*ctrl] );
}

int fstr_ctrl_seek_next_header_( int* ctrl ) {
	return c_fstr_ctrl_seek_next_header( ctrl_list[*ctrl] );
}

int fstr_ctrl_seek_next_header__( int* ctrl ) {
	return c_fstr_ctrl_seek_next_header( ctrl_list[*ctrl] );
}

int FSTR_CTRL_SEEK_NEXT_HEADER( int* ctrl ) {
	return c_fstr_ctrl_seek_next_header( ctrl_list[*ctrl] );
}

int FSTR_CTRL_SEEK_NEXT_HEADER_( int* ctrl ) {
	return c_fstr_ctrl_seek_next_header( ctrl_list[*ctrl] );
}

int FSTR_CTRL_SEEK_NEXT_HEADER__( int* ctrl ) {
	return c_fstr_ctrl_seek_next_header( ctrl_list[*ctrl] );
}

/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_c_h_name( int* ctrl, char* header_name, int *buff_size )
{
	char c_buff[STR_SIZE];

	if( c_fstr_ctrl_get_current_header_name( ctrl_list[*ctrl], c_buff )) return -1;

	strcpy_c2f( header_name, *buff_size, c_buff);
	return 0;
}

int fstr_ctrl_get_c_h_name_( int* ctrl, char* header_name, int *buff_size ) {
	return fstr_ctrl_get_c_h_name( ctrl, header_name, buff_size );
}

int fstr_ctrl_get_c_h_name__( int* ctrl, char* header_name, int *buff_size ) {
	return fstr_ctrl_get_c_h_name( ctrl, header_name, buff_size );
}

int FSTR_CTRL_GET_C_H_NAME( int* ctrl, char* header_name, int *buff_size ) {
	return fstr_ctrl_get_c_h_name( ctrl, header_name, buff_size );
}

int FSTR_CTRL_GET_C_H_NAME_( int* ctrl, char* header_name, int *buff_size ) {
	return fstr_ctrl_get_c_h_name( ctrl, header_name, buff_size );
}

int FSTR_CTRL_GET_C_H_NAME__( int* ctrl, char* header_name, int *buff_size ) {
	return fstr_ctrl_get_c_h_name( ctrl, header_name, buff_size );
}

/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_c_h_line_no( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_line_no( ctrl_list[*ctrl] );
}

int fstr_ctrl_get_c_h_line_no_( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_line_no( ctrl_list[*ctrl] );
}

int fstr_ctrl_get_c_h_line_no__( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_line_no( ctrl_list[*ctrl] );
}

int FSTR_CTRL_GET_C_H_LINE_NO( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_line_no( ctrl_list[*ctrl] );
}

int FSTR_CTRL_GET_C_H_LINE_NO_( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_line_no( ctrl_list[*ctrl] );
}

int FSTR_CTRL_GET_C_H_LINE_NO__( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_line_no( ctrl_list[*ctrl] );
}


/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_c_h_pos( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_pos( ctrl_list[*ctrl] );
}

int fstr_ctrl_get_c_h_pos_( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_pos( ctrl_list[*ctrl] );
}

int fstr_ctrl_get_c_h_pos__( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_pos( ctrl_list[*ctrl] );
}

int FSTR_CTRL_GET_C_H_POS( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_pos( ctrl_list[*ctrl] );
}

int FSTR_CTRL_GET_C_H_POS_( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_pos( ctrl_list[*ctrl] );
}

int FSTR_CTRL_GET_C_H_POS__( int* ctrl ) {
	return c_fstr_ctrl_get_current_header_pos( ctrl_list[*ctrl] );
}

/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_param( int* ctrl, const char* param_name, const char* value_list, char *type, void* val )
{
	int rcode;
	char p_name[STR_SIZE];
	char v_list[STR_SIZE];

	strcpy_f2c( p_name, param_name, STR_SIZE );
	strcpy_f2c( v_list, value_list, STR_SIZE );

	rcode = c_fstr_ctrl_get_param( ctrl_list[*ctrl], p_name, v_list, *type, val );

	if( rcode == 0 && (*type == 'S' || *type == 's') ){
		char buff[HECMW_NAME_LEN+1];
		strcpy( buff, (char*)val);
		strcpy_c2f( (char*)val, HECMW_NAME_LEN, buff );
	}
	return rcode;
}

int fstr_ctrl_get_param_( int* ctrl, const char* param_name, const char* value_list, char *type, void* val ) {
	return fstr_ctrl_get_param( ctrl, param_name, value_list, type, val );
}

int fstr_ctrl_get_param__( int* ctrl, const char* param_name, const char* value_list, char *type, void* val ) {
	return fstr_ctrl_get_param( ctrl, param_name, value_list, type, val );
}

int FSTR_CTRL_GET_PARAM( int* ctrl, const char* param_name, const char* value_list, char *type, void* val ) {
	return fstr_ctrl_get_param( ctrl, param_name, value_list, type, val );
}

int FSTR_CTRL_GET_PARAM_( int* ctrl, const char* param_name, const char* value_list, char *type, void* val ) {
	return fstr_ctrl_get_param( ctrl, param_name, value_list, type, val );
}

int FSTR_CTRL_GET_PARAM__( int* ctrl, const char* param_name, const char* value_list, char *type, void* val ) {
	return fstr_ctrl_get_param( ctrl, param_name, value_list, type, val );
}

/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_param_ex( int* ctrl,
		const char* param_name, const char* value_list, int *necessity, char *type, void* val )
{
	int rcode;
	char p_name[STR_SIZE];
	char v_list[STR_SIZE];

	strcpy_f2c( p_name, param_name, STR_SIZE );
	strcpy_f2c( v_list, value_list, STR_SIZE );

	rcode = c_fstr_ctrl_get_param_ex( ctrl_list[*ctrl], p_name, v_list, *necessity, *type, val );

	if( rcode_of_get_param == FSTR_CTRL_RCODE_PARAM_SUCCESS && (*type == 'S' || *type == 's') ){
		char buff[HECMW_NAME_LEN+1];
		strcpy( buff, (char*)val);
		strcpy_c2f( (char*)val, HECMW_NAME_LEN, buff );
	}
	return rcode;
}

int fstr_ctrl_get_param_ex_( int* ctrl,
		const char* param_name, const char* value_list, int *necessity, char *type, void* val )
{
	return fstr_ctrl_get_param_ex( ctrl, param_name, value_list, necessity, type, val );
}

int fstr_ctrl_get_param_ex__( int* ctrl,
		const char* param_name, const char* value_list, int *necessity, char *type, void* val )
{
	return fstr_ctrl_get_param_ex( ctrl, param_name, value_list, necessity, type, val );
}

int FSTR_CTRL_GET_PARAM_EX( int* ctrl,
		const char* param_name, const char* value_list, int *necessity, char *type, void* val )
{
	return fstr_ctrl_get_param_ex( ctrl, param_name, value_list, necessity, type, val );
}

int FSTR_CTRL_GET_PARAM_EX_( int* ctrl,
		const char* param_name, const char* value_list, int *necessity, char *type, void* val )
{
	return fstr_ctrl_get_param_ex( ctrl, param_name, value_list, necessity, type, val );
}

int FSTR_CTRL_GET_PARAM_EX__( int* ctrl,
		const char* param_name, const char* value_list, int *necessity, char *type, void* val )
{
	return fstr_ctrl_get_param_ex( ctrl, param_name, value_list, necessity, type, val );
}

/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_data_line_n( int* ctrl )
{
	return c_fstr_ctrl_get_data_line_n( ctrl_list[*ctrl] );
}

int fstr_ctrl_get_data_line_n_( int* ctrl )
{
	return fstr_ctrl_get_data_line_n( ctrl );
}

int fstr_ctrl_get_data_line_n__( int* ctrl ) {
	return fstr_ctrl_get_data_line_n( ctrl );
}

int FSTR_CTRL_GET_DATA_LINE_N( int* ctrl ) {
	return fstr_ctrl_get_data_line_n( ctrl );
}

int FSTR_CTRL_GET_DATA_LINE_N_( int* ctrl ) {
	return fstr_ctrl_get_data_line_n( ctrl );
}

int FSTR_CTRL_GET_DATA_LINE_N__( int* ctrl ) {
	return fstr_ctrl_get_data_line_n( ctrl );
}

/*-----------------------------------------------------------------------------------*/


int fstr_ctrl_get_data_n_in_line( int* ctrl, int *line_no, const char* delim )
{
	char delim_c[STR_SIZE];
	strcpy_f2c( delim_c, delim, STR_SIZE);

	return c_fstr_ctrl_get_data_n_in_line( ctrl_list[*ctrl], *line_no, delim_c );
}

int fstr_ctrl_get_data_n_in_line_( int* ctrl, int *line_no, const char* delim ) {
	return fstr_ctrl_get_data_n_in_line( ctrl, line_no, delim );
}

int fstr_ctrl_get_data_n_in_line__( int* ctrl, int *line_no, const char* delim ) {
	return fstr_ctrl_get_data_n_in_line( ctrl, line_no, delim );
}

int FSTR_CTRL_GET_DATA_N_IN_LINE( int* ctrl, int *line_no, const char* delim ) {
	return fstr_ctrl_get_data_n_in_line( ctrl, line_no, delim );
}

int FSTR_CTRL_GET_DATA_N_IN_LINE_( int* ctrl, int *line_no, const char* delim ) {
	return fstr_ctrl_get_data_n_in_line( ctrl, line_no, delim );
}

int FSTR_CTRL_GET_DATA_N_IN_LINE__( int* ctrl, int *line_no, const char* delim ) {
	return fstr_ctrl_get_data_n_in_line( ctrl, line_no, delim );
}


/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_data_error_pos(void) { return c_fstr_ctrl_get_data_error_pos(); }
int fstr_ctrl_get_data_error_pos_(void) { return c_fstr_ctrl_get_data_error_pos(); }
int fstr_ctrl_get_data_error_pos__(void) { return c_fstr_ctrl_get_data_error_pos(); }
int FSTR_CTRL_GET_DATA_ERROR_(void) { return c_fstr_ctrl_get_data_error_pos(); }
int FSTR_CTRL_GET_DATA_ERROR__(void) { return c_fstr_ctrl_get_data_error_pos(); }

int fstr_ctrl_get_data_error_line(void) { return c_fstr_ctrl_get_data_error_line(); }
int fstr_ctrl_get_data_error_line_(void) { return c_fstr_ctrl_get_data_error_line(); }
int fstr_ctrl_get_data_error_line__(void) { return c_fstr_ctrl_get_data_error_line(); }
int FSTR_CTRL_GET_DATA_ERROR_LINE(void) { return c_fstr_ctrl_get_data_error_line(); }
int FSTR_CTRL_GET_DATA_ERROR_LINE_(void) { return c_fstr_ctrl_get_data_error_line(); }
int FSTR_CTRL_GET_DATA_ERROR_LINE__(void) { return c_fstr_ctrl_get_data_error_line(); }

/*-----------------------------------------------------------------------------------*/



int fstr_ctrl_get_data_v_f( fstr_ctrl_data* ctrl, int line_no, const char* format, va_list va )
{
	char fmt_c[ STR_SIZE ];

	strcpy_f2c( fmt_c, format, STR_SIZE);
	return c_fstr_ctrl_get_data_v( ctrl, line_no, format, va );
}


int fstr_ctrl_get_data( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;

	va_start(va, format);
	r = fstr_ctrl_get_data_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);

	return r;
}

int fstr_ctrl_get_data_( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;

	va_start(va, format);
	r = fstr_ctrl_get_data_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);

	return r;
}

int fstr_ctrl_get_data__( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;

	va_start(va, format);
	r = fstr_ctrl_get_data_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);

	return r;
}

int FSTR_CTRL_GET_DATA( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;

	va_start(va, format);
	r = fstr_ctrl_get_data_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);

	return r;
}

int FSTR_CTRL_GET_DATA_( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;

	va_start(va, format);
	r = fstr_ctrl_get_data_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);

	return r;
}

int FSTR_CTRL_GET_DATA__( int* ctrl, int* line_no, const char* format, ... )
{
	va_list va;
	int r;

	va_start(va, format);
	r = fstr_ctrl_get_data_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);

	return r;
}
/*-----------------------------------------------------------------------------------*/



int fstr_ctrl_get_data_ex_v_f( fstr_ctrl_data* ctrl, int line_no, const char* format, va_list va )
{
	int r;
	char fmt_c[ STR_SIZE ];

	strcpy_f2c( fmt_c, format, STR_SIZE);
	r = c_fstr_ctrl_get_data_v( ctrl, line_no, fmt_c, va );

	if( r != 0 ) {
		set_record_data_line_err_msg( ctrl, r);
		return -1;
	}

	return 0;
}

int fstr_ctrl_get_data_ex( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_ex_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);
	return r;
}

int fstr_ctrl_get_data_ex_( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_ex_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);
	return r;
}

int fstr_ctrl_get_data_ex__( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_ex_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);
	return r;
}

int FSTR_CTRL_GET_DATA_EX( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_ex_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);
	return r;
}

int FSTR_CTRL_GET_DATA_EX_( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_ex_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);
	return r;
}

int FSTR_CTRL_GET_DATA_EX__( int* ctrl, int *line_no, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_ex_v_f( ctrl_list[*ctrl], *line_no, format, va );
	va_end(va);
	return r;
}



/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_get_data_array_ex_v_f( fstr_ctrl_data* ctrl, const char* format, va_list va )
{
	int r;
	char fmt_c[ STR_SIZE ];

	strcpy_f2c( fmt_c, format, STR_SIZE);

	fg_fortran_get_data_array_v = 1;
	r = c_fstr_ctrl_get_data_array_v( ctrl, fmt_c, va );
	fg_fortran_get_data_array_v = 0;

	if( r != 0 ) {
		set_record_data_line_err_msg( ctrl, r);
		return -1;
	}

	return 0;
}

int fstr_ctrl_get_data_array_ex( int* ctrl, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_array_ex_v_f( ctrl_list[*ctrl], format, va );
	va_end(va);
	return r;
}

int fstr_ctrl_get_data_array_ex_( int* ctrl, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_array_ex_v_f( ctrl_list[*ctrl], format, va );
	va_end(va);
	return r;
}

int fstr_ctrl_get_data_array_ex__( int* ctrl, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_array_ex_v_f( ctrl_list[*ctrl], format, va );
	va_end(va);
	return r;
}

int FSTR_CTRL_GET_DATA_ARRAY_EX( int* ctrl, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_array_ex_v_f( ctrl_list[*ctrl], format, va );
	va_end(va);
	return r;
}

int FSTR_CTRL_GET_DATA_ARRAY_EX_( int* ctrl, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_array_ex_v_f( ctrl_list[*ctrl], format, va );
	va_end(va);
	return r;
}

int FSTR_CTRL_GET_DATA_ARRAY_EX__( int* ctrl, const char* format, ... )
{
	va_list va;
	int r;
	va_start(va, format);
	r = fstr_ctrl_get_data_array_ex_v_f( ctrl_list[*ctrl], format, va );
	va_end(va);
	return r;
}


/*-----------------------------------------------------------------------------------*/

int fstr_ctrl_close( int* ctrl )
{
	int fg = c_fstr_ctrl_close( ctrl_list[*ctrl] );
	HECMW_free(ctrl_list[*ctrl]);
	ctrl_list[*ctrl] = NULL;
	return fg;
}


int fstr_ctrl_close_( int* ctrl ) { return fstr_ctrl_close( ctrl ); }
int fstr_ctrl_close__( int* ctrl ) { return fstr_ctrl_close( ctrl ); }
int FSTR_CTRL_CLOSE( int* ctrl )  { return fstr_ctrl_close( ctrl ); }
int FSTR_CTRL_CLOSE_( int* ctrl )  { return fstr_ctrl_close( ctrl ); }
int FSTR_CTRL_CLOSE__( int* ctrl ) { return fstr_ctrl_close( ctrl ); }

/*-----------------------------------------------------------------------------------*/

void fstr_ctrl_dump(  int* ctrl ) {
	c_fstr_ctrl_dump(  ctrl_list[*ctrl] );
}

void fstr_ctrl_dump_(  int* ctrl ) {
	c_fstr_ctrl_dump(  ctrl_list[*ctrl] );
}

void fstr_ctrl_dump__(  int* ctrl ) {
	c_fstr_ctrl_dump(  ctrl_list[*ctrl] );
}

void FSTR_CTRL_DUMP(  int* ctrl ) {
	c_fstr_ctrl_dump(  ctrl_list[*ctrl] );
}

void FSTR_CTRL_DUMP_(  int* ctrl ) {
	c_fstr_ctrl_dump(  ctrl_list[*ctrl] );
}

void FSTR_CTRL_DUMP__(  int* ctrl ) {
	c_fstr_ctrl_dump(  ctrl_list[*ctrl] );
}

/* ==================================================================================*/
/* LOCAL UTIRITY FUNCTION                                                            */
/* ==================================================================================*/


/* JP-20 */

static
void set_param_err_msg( fstr_ctrl_data* ctrl, const char* param_name, const char* msg )
{
	sprintf(err_msg, "fstr control file error(param): line:%d, %s, %s\n",
		c_fstr_ctrl_get_current_header_line_no( ctrl ),
		param_name, msg );
}


/* JP-21 */

static
void set_data_err_msg( fstr_ctrl_data* ctrl, int line, int data_no, const char* msg )
{
	line += c_fstr_ctrl_get_current_header_line_no( ctrl );
	sprintf( err_msg, "fstr control file error(data): line:%d, column:%d : %s\n", line, data_no, msg);
}



static void set_record_data_line_err_msg(  fstr_ctrl_data* ctrl, int r )
{
	char msg[buffsize];

	int line_no = fstr_ctrl_get_data_error_line();
	int pos = fstr_ctrl_get_data_error_pos();

	switch( r ) {
	case 0:
		strcpy( msg, "no error");
		break;
	case FSTR_CTRL_RCODE_DATA_TYPE_ERROR:
		strcpy( msg, "data type converting error");
		break;
	case FSTR_CTRL_RCODE_DATA_RANGE_ERROR:
		strcpy( msg, "data range error");
		break;
	case FSTR_CTRL_RCODE_DATA_NOTHING:
		strcpy( msg,"data must exist" );
		break;
	case FSTR_CTRL_RCODE_DATA_LINE_NOTHING:
		strcpy( msg,"data line does not exist" );
		break;
	default:
		sprintf( msg, "data line unknown error (r:%d)", r);
	}
	set_data_err_msg( ctrl, line_no, pos, msg );
}


/* JP-22 */

static
char* gettoken( const char* line )
{
	static char buff[buffsize];
	static char* p;
	static char* h;
	char* token;
	char* t;
	int is_null = 0;

	if( line ) {
		strcpy( buff, line );
		h = p = buff;
	} else {
		if( *p == 0 || *p == '\n' || *p == '\r'  )
			return NULL;
	}
	for(;;){
		if( *p == 0 || *p == '\n' || *p == '\r'  || *p == ',' ) {
			if( *p == 0 ){
				is_null = 1;
			} 
			*p = 0;
			/* remove space in head */
			while( *h == ' ' ) h++;
			/* remove space in tail */
			t = p;
			if( t != buff){
				t--;
				while( *t == ' ' && t != h ) {
					*t = 0;
					t--;
				}
			}
			token = h;
			h = p;
			h++;
			break;
		}
		p++;
	}
	if( !is_null )
		p++;
	return token;
}


/* JP-23 */
/* JP-24 */
/* JP-25 */
/* JP-26 */

static
void strcpy_f2c( char* dest, const char* src, int len)
{
	int i;
	int fg = 0;

	for(i=0; i<len-1; i++){
		if( src[i] != ' ')
			fg = 1;
		else if( fg ){
			break;
		}
		dest[i] = src[i];
	}
	dest[i] = 0;
}



/* JP-27 */

static
void strcpy_c2f( char* dest, int len, const char* src)
{
	int i;
	int fg = 0;
	for(i=0; i<len; i++){
		if( src[i] == 0) fg = 1;
		if(fg)
			dest[i] = ' ';
		else
			dest[i] = src[i];
	}
}

/* JP-28 */

static char* remove_header_space( char* token )
{
	char* p = token;
	while( *p ){
		if( *p != ' ') break;
		p++;
	}
	return p;
}


/* JP-29 */

static
int Strncmpi( const char* s1, const char* s2, int len)
{
	char *p1, *p2;
	int i;

	p1 = (char*)s1;
	p2 = (char*)s2;

	i = 0;
	while( *p1 && *p2 && i<len ){
		if(toupper( *p1 ) != toupper( *p2 ))
			return -1;
		p1++;
		p2++;
		i++;
	}
	return 0;
}

static
int Strcmpi( const char* s1, const char* s2 )
{
	char *p1, *p2;

	p1 = (char*)s1;
	p2 = (char*)s2;

	while( *p1 && *p2 ){
		if(toupper( *p1 ) != toupper( *p2 ))
			return -1;
		p1++;
		p2++;
	}
	if( *p1 || *p2 )
		return -1;
	else
		return 0;
}

/* JP-30 */
#if 0
static
void Strupr( char* s )
{
	while( *s ){
		*s = toupper(*s);
		s++;
	}
}
#endif

/* JP-31 */

static void remove_cr( char* s )
{
	while( *s ) {
		if( *s == '\n' || *s == '\r' ) {
			*s = 0;
			return;
		}
		s++;
	}
}



/* JP-32 */

static void format_conv( const char* format, char* fmt, int* array_size )
{
	char* p = (char*)format;
	int i = 0;
	int num_fg = 0;
	char num_str[buffsize];
	int num_index = 0;
	int n;

	while(*p){
		if( '0' <= *p && *p <= '9' ) {
			num_fg = 1;
			num_str[ num_index ] = *p;
			num_index++;
		} else {
			if( num_fg ){
				if( i == 0 ) assert(0);
				num_str[num_index] = 0;
				sscanf( num_str, "%d", &n);
				array_size[i-1] = n;
				num_fg = 0;
				num_index = 0;
			}
			fmt[i] = *p;
			array_size[i] = 1;
			i++;
		}
		p++;
	}
	fmt[i] = 0;
	if( num_fg ){
		if( i == 0 ) assert(0);
		num_str[num_index] = 0;
		sscanf( num_str, "%d", &n);
		array_size[i-1] = n;
	}
}
