/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Noboru Imai (Univ. of Tokyo)                  *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

/* JP-0 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_io_mesh.h"
#include "hecmw_io_struct.h"
#include "hecmw_struct.h"
#include "hecmw_config.h"
#include "hecmw_system.h"
#include "hecmw_dist.h"
#include "hecmw_dist_print.h"
#include "hecmw_common.h"
#include "hecmw_path.h"
#include "hecmw_conn_conv.h"
#include "hecmw_io_nastran.h"

/* def.h */ 


/* JP-1 */


#ifndef NULL
#define NULL ((void*)0)
#endif


/* JP-2 */


#define FIELD_BUFFER_SIZE  16+1 			/* JP-3 */
#define FIELD_NUMBER  10				/* JP-4 */
#define LINE_BUFFER_SIZE  (FIELD_BUFFER_SIZE * FIELD_NUMBER) /* JP-5 */

#define FIRST_FIELD_COLUMN 8	/* JP-6 */
#define LAST_FIELD_COLUMN 8	/* JP-7 */
#define SMALL_FIELD_COLUMN 8	/* JP-8 */
#define LARGE_FIELD_COLUMN 16	/* JP-9 */
#define LARGE_FIELD_LINE_FIELD_NUMBER 6	/* JP-10 */

/* JP-11 */


#define IFF_FIELD_NUMBER FIELD_NUMBER
#define IFF_FIELD_SIZE 16				/* JP-12 */
#define IFF_FIELD_BUFFER_SIZE  (IFF_FIELD_SIZE+1) 	/* JP-13 */
#define IFF_SIZE (IFF_FIELD_BUFFER_SIZE * IFF_FIELD_NUMBER)	/* JP-14 */

#define IFF_LAST_FIELD 9	/* JP-15 */


/* JP-16 */

typedef struct ST_IFF_NODE {
	char iff[IFF_SIZE];
	int line_no;			/* JP-17 */
	struct ST_IFF_NODE *next;
} iff_node_t;


/* JP-18 */

typedef struct ST_IFF_BULK {
	iff_node_t* iff_node_root;	/* JP-19 */
	iff_node_t* iff_node_last;
	int file_no;			/* JP-20 */
	struct ST_IFF_BULK* next_bulk;
} iff_bulk_t;


/* JP-21 */

typedef struct ST_IFF_BULK_LIST {
	iff_bulk_t* bulk_root;
	int bulk_number;
} iff_bulk_list_t;


/* JP-22 */

typedef enum ENUM_FIELD_FORMAT {
	fft_none,		/* JP-23 */
	fft_comment,		/* JP-24 */
	fft_free_field_format,	/* JP-25 */
	fft_small_field_format,	/* JP-26 */
	fft_large_field_format	/* JP-27 */
} field_format_t;


/* JP-28 */


#ifndef HECMW_FILENAME_LEN
#define HECMW_FILENAME_LEN 1023
#endif


#define FILE_STACK_SIZE 10

typedef struct ST_FILE_STACK {
	FILE* fp;
	char filename[HECMW_FILENAME_LEN+1];
	int lineno;
} file_stack_t;


static file_stack_t file_stack[ FILE_STACK_SIZE ];
static int file_stack_pos;

static char grid_filename[HECMW_FILENAME_LEN+1] = "Unknown";

static FILE* c_fp = NULL;
static char c_filename[HECMW_FILENAME_LEN+1] = "Unknown";
static int c_lineno = 0;


/* JP-29 */


/* utils.c */

/* JP-30 */
static char* str_upr( char* s );

/* JP-31 */
static
void set_error_of_field_convert( const char* file_name, int line_no, int field_no,
		const char* bulk_name, const char* param_name, const char* param_val, int type );

/* JP-32 */
static
void set_error_of_blank_field( const char* file_name, int line_no, int field_no,
				const char* bulk_name, const char* param_name);

/* JP-33 */
static void remove_cr( char* line );

/* JP-34 */
static void skip_to_begin_bulk( FILE* fp, int* line_no);


/* JP-35 */
static char *ngrp_name_by_GID(int gid, char* name);

/* JP-36 */
static char *egrp_name_by_PID(int pid, char* name);

/* JP-37 */
static char *egrp_name_by_SID(int sid, char* name);

/* JP-38 */
static char *ngrp_name_by_SID_GID( int sid, int gid, char* name);

/* JP-39 */
static char *ngrp_name_by_SID(int sid, int sub_id, char* name);

/* JP-40 */
static char *matrial_name_by_MID(int mid, char* name);

/* JP-41 */
static char *egrp_CTRIAX6_name_by_MID(int mid, char* name);

/* JP-42 */
static int is_CTRIAX6_egrp_name( char* name, int* mid );


/* ------------------------------------------------------------------------------------------------ */
/* field_format_judge.c */

/* JP-43 */
static field_format_t field_format_judge( char* line );


/* ------------------------------------------------------------------------------------------------ */
/* iff.c */


static void iff_init( char* iff );				/* JP-44 */
static void iff_clear( char* iff );				/* JP-45 */
static int iff_field_empty( char* field );			/* JP-46 */
static int iff_is_blank_field( char* iff, int no );		/* JP-47 */
static char* iff_get_field( char* iff, int no );		/* JP-48 */
static void iff_set_field( char* iff, int no, char* field );	/* JP-49 */
static void iff_regulize( char* iff);				/* JP-50 */
static void iff_copy( char* dest, char* src );			/* JP-51 */
static char* iff_get_continuous_line_pointer( char* iff);	/* JP-52 */
static char* iff_get_pointing_line_pointer( char* iff );	/* JP-53 */

static int iff_is_parent( char* iff );				/* JP-54 */
static int iff_is_last( char* iff );				/* JP-55 */
static int iff_is_continuous_line( char* iff);			/* JP-56 */
static int iff_is_pointing_line( char* iff);			/* JP-57 */

static int iff_add_auto_pointer( char* iff1, char* iff2, int* counter );/* JP-58 */

/* JP-59 */
/* JP-60 */

static int iff_field_to_int( char* field, int default_val, int* val );				/* JP-61 */
static int iff_field_to_uint( char* field, unsigned int default_val, unsigned int* val );	/* JP-62 */
static int iff_field_to_double( char* field, double default_val, double* val );			/* JP-63 */
static int iff_field_to_param( char* field, char type, void* param, void** next_param_pos );	/* JP-64 */

/* JP-65 */

static void iff_dump(char* iff);
static void iff_dump_csv(char* iff);


/* ------------------------------------------------------------------------------------------------ */
/* iff_conv.c */

static void free_to_iff_format( char* line, char* iff);			/* JP-66 */
static void small_to_iff_format( char* line, char* iff);		/* JP-67 */
static void large_to_iff_format( char* line1, char* line2, char* iff);	/* JP-68 */

static char* get_fixed_token( char* src, char* token, int field_size ); /* JP-69 */


/* ------------------------------------------------------------------------------------------------ */
/* iff_node.c */

static iff_node_t* iff_node_create(void);				/* JP-70 */
static void iff_node_free( iff_node_t* node );				/* JP-71 */
static void iff_node_init( iff_node_t* node );				/* JP-72 */
static void iff_node_set( iff_node_t* node, char* iff, int line_no);	/* JP-73 */


/* ------------------------------------------------------------------------------------------------ */
/* read_iff.c */

/* JP-74 */
static int read_iff( FILE* fp, char* iff, int* line_no );

/* ------------------------------------------------------------------------------------------------ */
/* bulk.c */

static void iff_bulk_init( iff_bulk_t* bulk );	/* JP-75 */

/* JP-76 */

static iff_bulk_t* iff_bulk_create(void);						/* JP-77 */
static void iff_bulk_free( iff_bulk_t* bulk );					/* JP-78 */
static void iff_bulk_regist( iff_bulk_t* bulk, char* iff, int line_no );	/* JP-79 */
static void iff_bulk_move_append( iff_bulk_t* dest, iff_bulk_t* src );		/* JP-80 */
static int iff_bulk_is_completed( iff_bulk_t* bulk );				/* JP-81 */
static int iff_bulk_parse( iff_bulk_t*  bulk );					/* JP-82 */

/* JP-83 */

static char* iff_bulk_get_name( iff_bulk_t* bulk );				/* JP-84 */
static char* iff_bulk_get_field( iff_bulk_t* bulk, int field_no );		/* JP-85 */
static char* iif_bulk_get_param( iff_bulk_t* bulk, int param_no );		/* JP-86 */
static char* iff_bulk_get_param_data( iff_bulk_t* bulk, int param_no, int* line_no, int* field_no ); /* JP-87 */

/* JP-88 */
static int iff_bulk_get_param_list( iff_bulk_t* bulk, const char* format, int* result, ... );

/* JP-89 */
static int iff_bulk_get_param_list_pattern( iff_bulk_t* bulk, const char* format, int start_param_no, int* size, ... );


static int iff_bulk_get_line_number( iff_bulk_t* bulk );			/* JP-90 */
static int iff_bulk_get_line_no( iff_bulk_t* bulk, int iff_pos );		/* JP-91 */

/* JP-92 */

static iff_node_t* iff_bulk_search_pointing_line( iff_bulk_t* bulk, char* line_pointer );	/* JP-93 */
static iff_node_t* iff_bulk_search_continuous_line( iff_bulk_t* bulk, char* line_pointer );	/* JP-94 */

/* JP-95 */

static void iff_bulk_dump( iff_bulk_t* bulk);	/* JP-96 */

/* ------------------------------------------------------------------------------------------------ */
/* bulk_parse_head.c */

static char* iff_bulk_get_name( iff_bulk_t* bulk );			/* JP-97 */
static char* iff_bulk_get_param( iff_bulk_t* bulk, int field_no );	/* JP-98 */
static int iff_bulk_type( const char* bulkname );	/* JP-99 */

/* JP-100 */
static int iff_bulk_parse( iff_bulk_t* bulk );

/* ------------------------------------------------------------------------------------------------ */
/* bulk_parse.c */

/* JP-101 */
static struct hecmw_io_material * create_mat_struct(const char* name, int item_n, int sub_n[], double data[] );

/* JP-102 */
static int surface_elem_store( int bulk_type, unsigned int EID, unsigned int PID,
	unsigned int G[], int G_number, double THETA_MCID, double ZOFFS, double T[], int T_number);

/* JP-103 */
static int iff_bulk_parse_solid_elem( iff_bulk_t* bulk,
				int first_etype, int second_etype,
				int conv_table[],
				int g_number, int g_ness_number );

/* ------------------------------------------------------------------------------------------------ */
/* bulk_list.c */

static void iff_bulk_list_init( iff_bulk_list_t* list);	/* JP-104 */
static void iff_bulk_list_clear( iff_bulk_list_t* list);/* JP-105 */

static int iff_bulk_list_is_empty( iff_bulk_list_t* list); /* JP-106 */

/* JP-107 */
static void iff_bulk_list_regist( iff_bulk_list_t* bulk_list, iff_bulk_t* bulk);/* JP-108 */
static iff_bulk_t* iff_bulk_list_search_pointing_bulk( iff_bulk_list_t* bulk_list, char* line_pointer); /* JP-109 */
static iff_bulk_t* iff_bulk_list_search_continuous_bulk( iff_bulk_list_t* bulk_list, char* line_pointer);/* JP-110 */
static int iff_bulk_list_remove( iff_bulk_list_t* bulk_list, iff_bulk_t* bulk );/* JP-111 */

/* JP-112 */
static int iff_bulk_list_after_operation( iff_bulk_list_t* list );

/* JP-113 */
static void iff_bulk_list_dump( iff_bulk_list_t * bulk_list ); /* JP-114 */

/* ------------------------------------------------------------------------------------------------ */
/* iff_operation.c */

/* JP-115 */
static int iff_operation( iff_bulk_list_t* bulk_list, char* iff, int line_no);

/* ------------------------------------------------------------------------------------------------ */
/* f_open_close.c */

/* JP-116 */
static void file_stack_init(void);				/* JP-117 */
static int file_stack_push(char* fname, FILE* fp, int lineno);	/* JP-118 */
static int file_stack_pop(char* fname, FILE** fp, int* lineno);	/* JP-119 */
static int file_stack_clear(void); 				/* JP-120 */
static int file_stack_check( const char* fname );		/* JP-121 */

/* JP-122 */
static void nastran_file_init(void);			/* JP-123 */
static int nastran_file_open( const char* filename );	/* JP-124 */
static int nastran_file_close(void);			/* JP-125 */

/* ------------------------------------------------------------------------------------------------ */
/* read_nastran.c */

/* JP-126 */
static int read_nastran( const char* filename );

/* list_creator.c */ 

/* JP-127 */



/* JP-128 */


#define CREATE_LIST_CODE( type ) \
\
static int is_eq_##type ( type a, type b);\
static void copy_##type ( type * dest, type src );\
static void delete_##type ( type a );\
\
typedef struct struct_list_node_##type { \
	type data; \
	struct struct_list_node_##type *next; \
} t_list_node_##type ; \
\
\
typedef struct struct_list_##type { \
	t_list_node_##type *root;\
	t_list_node_##type *current, *next;\
} t_list_##type ;\
\
\
static void list_init_##type ( t_list_##type *list ){\
	list->root = NULL; \
} \
\
static void list_clear_##type ( t_list_##type *list ){\
	t_list_node_##type *t;\
	t_list_node_##type *node = list->root;\
	while( node ) {\
		t = node->next;\
		delete_##type ( node->data );\
		node = t;\
	}\
	list->root = NULL; \
} \
\
\
static t_list_node_##type *list_add_##type ( t_list_##type *list, type d ){\
	t_list_node_##type *node;\
	node = HECMW_malloc( sizeof( t_list_node_##type ));\
	copy_##type( &node->data, d); \
	node->next = NULL; \
	if( list->root ) { \
		node->next = list->root; \
		list->root = node; \
	} else { \
		list->root = node; \
	} \
	return node; \
}\
\
\
\
static t_list_node_##type *list_search_##type ( t_list_##type *list, type data ) {\
	static t_list_##type *list_s; \
	if( list != NULL) {\
		list->current = list->root; \
		list_s = list; \
	} \
	while( list_s->current ) { \
		if( is_eq_##type( list_s->current->data, data)) \
			return list_s->current; \
		list_s->current = list_s->current->next; \
	}\
	return NULL; \
}\
\
\
static void list_scan_start_##type ( t_list_##type *list ) {\
	list->next = list->current = list->root; \
}\
\
\
static int list_scan_next_##type ( t_list_##type *list ) {\
	t_list_node_##type *t;\
	if( list->next == NULL )\
		return 0;\
\
	t = list->next;\
	list->next = list->next->next;\
	list->current = t;\
\
	return 1;\
}\
\
\
static void list_get_current_data_##type ( t_list_##type *list, type *data ) {\
	*data = list->current->data;\
}\
\
\
static t_list_node_##type * list_get_current_node_##type ( t_list_##type *list ) {\
	return list->current; \
}\
\
\
static int list_operation_##type ( t_list_##type *list, int (*foo)( type )) {\
	list_scan_start_##type ( (list) );\
	while( list_scan_next_##type ((list)) ) {\
		int err =  (*(foo))( list->current->data );\
		if( err ) return err;\
	}\
	return 0;\
}

/* utils.c */ 

/* JP-129 */


/* JP-130 */

static
char* str_upr( char* s )
{
	char* p = s;
	while( *p ){
		if( 'a' <= *p  && *p <'A' )
			*p += ('A' - 'a');
		p++;
	}
	return s;
}

/* JP-131 */

/* JP-132 */

static
void set_error_of_field_convert( const char* file_name, int line_no, int field_no,
				const char* bulk_name, const char* param_name, const char* param_val, int type )
{
	char line[256];
	char msg[256];
	char type_s[10];
	char val_str[256];

	switch(type){
	case 's': case'S':
		strcpy( type_s, "string");
		break;
	case 'i': case'I':
	case 'u': case'U':
		strcpy( type_s, "integer");
		break;
	case 'd': case'D':
		strcpy( type_s, "real");
		break;
	case '_':
		strcpy( type_s, "blank");
		break;
	default:
		HECMW_assert(0);
	}

	if( param_val ){
		int i;
		char s[20];
		strcpy(val_str, "H");
		for(i=0; i<IFF_FIELD_SIZE && param_val[i]; i++) {
			sprintf(s, "%x", param_val[i]);
			strcat( val_str, s);
		}
	}else
		strcpy( val_str, "NULL");

	sprintf(line, "%s:%d:", file_name, line_no );
	sprintf(msg, "Cannot convert '%s'(%s) to '%s' value ('%s' (field.%d) in '%s')",
				param_val, val_str, type_s, param_name, field_no, bulk_name );

	HECMW_set_error( HECMW_IO_HEC_E0001, "%s%s", line, msg );/* JP-133 */
}


/* JP-134 */

static
void set_error_of_blank_field( const char* file_name, int line_no, int field_no,
				const char* bulk_name, const char* param_name)
{
	char line[256];
	char msg[256];

	sprintf(line, "%s:%d:", file_name, line_no );
	sprintf(msg, "Cannot be blank field ('%s' (no.%d) in '%s')", param_name, field_no, bulk_name );

	HECMW_set_error( HECMW_IO_HEC_E0001, "%s%s", line, msg );/* JP-135 */
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-136 */

#define REMOVE_CR_MAX_LINE 256

static
void remove_cr( char* line )
{
	int i = 0;
	char* p = line;
	while( *p && i<REMOVE_CR_MAX_LINE ){
		if(*p == '\n' || *p == 0xd ){
			*p = 0;
			return;
		}
		i++;
		p++;
	}
}



/* ------------------------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------------------------ */
/* JP-137 */

static
void skip_to_begin_bulk( FILE* fp, int* line_no)
{
	char line[ IFF_SIZE ];
	
	while( fgets( line, IFF_SIZE, fp) ){
		(*line_no)++;
		if(strncmp( line, "BEGIN BULK", sizeof("BEGIN BULK")-1) == 0)
			return;
	}
}



/* JP-138 */



/* JP-139 */

static char *ngrp_name_by_GID(int gid, char* name)
{
	sprintf( name, "G%d", gid);
	return name;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-140 */

static char *egrp_name_by_PID(int pid, char* name)
{
	sprintf( name, "P%d", pid);
	return name;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-141 */

static char* create_SID_grp_name( char* name )
{
	static int SID_counter = 1;
	sprintf( name, "SID%d", SID_counter);
	SID_counter++;
	return name;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-142 */


static char *egrp_name_by_SID(int sid, char* name)
{
	/* sprintf( name, "S%d", sid); */
	create_SID_grp_name(name);
	return name;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-143 */

static char *ngrp_name_by_SID_GID( int sid, int gid, char* name)
{
	/* sprintf( name, "S%dG%d", sid, gid); */
	create_SID_grp_name(name);
	return name;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-144 */

static char *ngrp_name_by_SID(int sid, int sub_id, char* name)
{
	/* sprintf( name, "SID%d-%d", sid, sub_id); */
	create_SID_grp_name(name);
	return name;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-145 */

static char *matrial_name_by_MID(int mid, char* name)
{
	sprintf( name, "M%d", mid);
	return name;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-146 */

#define CTRIAX6_EGRP_NAME_HEADER "CTRIAX6-"

static char *egrp_CTRIAX6_name_by_MID(int mid, char* name)
{
	sprintf( name, "%s%d", CTRIAX6_EGRP_NAME_HEADER, mid);
	return name;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-147 */

static int is_CTRIAX6_egrp_name( char* name, int* mid )
{
	char* p = CTRIAX6_EGRP_NAME_HEADER;
	char* s = name;
	while( *s && *p ){
		if( *s != *p)
			return 0;
		s++;
		p++;
	}
	*mid = atoi(s);
	return 1;
}


/* field_format_judge.c */ 

/* JP-148 */

/* JP-149 */


/* JP-150 */

/* ------------------------------------------------------------------------------------------------ */
/* JP-151 */
/* JP-152 */

static
field_format_t field_format_judge( char* line )
{
	int i;
	int fg;
	int n;
	char* p;

	p = line;

	if( *p == '$' || *p == 0)
		return fft_comment;

	fg = 1;
	while( *p != 0){
		if(*p == ',')
			return fft_free_field_format;
		if(*p != ' ')
			fg = 0;
		p++;
	}

	if(fg)
		return fft_comment;

	n = strlen( line );
	p = line;

	if( n < FIRST_FIELD_COLUMN )
		p += n;
	else
		p += FIRST_FIELD_COLUMN;

	while( p!=line){
		p--;
		if( *p != ' '){
			if( *p == '*')
				return fft_large_field_format;
			else
				return fft_small_field_format;
		}
	}

	return 	fft_small_field_format;
}

/* iff.c */ 

/* JP-153 */



/* JP-154 */

static
void iff_init( char* iff )
{
	iff_clear( iff );
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-155 */

static
void iff_clear( char* iff )
{
	memset( iff, 0, IFF_SIZE );
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-156 */

static
int iff_field_empty( char* field )
{
	return (field[0] == 0);
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-157 */

static int iff_is_blank_field( char* iff, int no )
{
	char* field = iff_get_field(iff, no);
	return field[0] == 0;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-158 */

static
char* iff_get_field( char* iff, int no )
{
	char *p = iff;

	if( no<0 || no>= IFF_FIELD_NUMBER )
		return NULL;

	p += IFF_FIELD_BUFFER_SIZE * no;
	return p;
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-159 */
/* JP-160 */



static
void iff_set_field( char* iff, int no, char*  field )
{
	int i, fg;
	char *start_p, *p, *s;

	HECMW_assert( no < FIELD_NUMBER );

	s = field;

	start_p = iff;
	start_p += IFF_FIELD_BUFFER_SIZE * no;

	/* JP-161 */

	fg = 0;
	p = start_p;
	for( i=0; *s && i< IFF_FIELD_SIZE; i++, s++ ){
		if( *s != ' ')
			fg = 1;

		if( fg ){
			*p = *s;
			p++;
		}
	}
	*p = 0;

	/* JP-162 */

	if(p == start_p)
		return;

	do {
		p--;
		if(*p == ' ')
			*p = 0;
		else
			break;
	} while( p != start_p );
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-163 */


static void iff_regulize( char* iff)
{
	char *last_field;

	/* JP-164 */

	int i;
	char* p = iff;
	p += IFF_FIELD_SIZE;

	do {
		p--;
		if(*p == '*'){
			*p = 0;
			break;
		}
	} while( p != iff );

	/* JP-165 */

	if(p == iff)
		return;

	do {
		p--;
		if(*p == ' ')
			*p = 0;
		else
			break;
	} while( p != iff );

	/* JP-166 */

	/* JP-167 */
	/* JP-168 */

	if( iff[0] == '*' ){
		if( iff[1] == 0 )
			iff[0] = 0;
		else
			iff[0] = '+';
	}

	last_field = iff_get_field( iff, IFF_LAST_FIELD );
	if( last_field[0] == '*' ){
		if(last_field[1] == 0)
			last_field[0] = 0;
		else
			last_field[0] = '+';
	}
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-169 */

static
void iff_copy( char* dest, char* src )
{
	memcpy( dest, src, IFF_SIZE);
}


/* JP-170 */

/* ------------------------------------------------------------------------------------------------ */
/* JP-171 */

static
int iff_is_continuous_line( char* iff)
{
	char* f = iff_get_field( iff, 0 );

	return (f[0] == '+');
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-172 */

static
int iff_is_pointing_line( char* iff)
{
	char* f = iff_get_field( iff, IFF_LAST_FIELD );
	return !iff_field_empty(f);
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-173 */

static
char* iff_get_continuous_line_pointer( char* iff)
{
	char* f = iff_get_field( iff, 0 );

	if(iff_field_empty(f))
		return NULL;

	if(f[0] == '+'){
		if( f[1] == '+') /* JP-174 */
			f++;
		return f;
	} else
		return NULL;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-175 */
/* JP-176 */

static
char* iff_get_pointing_line_pointer( char* iff )
{
	char* f = iff_get_field( iff, IFF_LAST_FIELD );
	if(iff_field_empty(f))
		return NULL;
	else {
		return f;
	}
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-177 */

static
int iff_is_parent( char* iff )
{
	char* field = iff_get_field( iff, 0 );
	return (field[0] != '+');
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-178 */

static
int iff_is_last( char* iff )
{
	return (iff_get_pointing_line_pointer( iff ) == NULL );
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-179 */
/* JP-180 */


static
int iff_add_auto_pointer( char* iff1, char* iff2, int* counter )
{
	if( iff_is_blank_field( iff2, 0 )){
		if( iff_is_blank_field( iff1, IFF_LAST_FIELD )){
			char pointer[IFF_FIELD_SIZE];
			sprintf( pointer, "+%d", *counter);
			(*counter)++;
			iff_set_field( iff1, IFF_LAST_FIELD, pointer);
			iff_set_field( iff2, 0, pointer);
		}else
			return -1;
	}
	return 0;
}




/* ------------------------------------------------------------------------------------------------ */
/* JP-181 */

static
int iff_field_to_uint( char* field, unsigned int default_val, unsigned int* val )
{
	char* endptr;

	if( field == NULL || (field && *field == 0)) {
		*val = default_val;
		return 2;
	}

	*val = strtoul( field, &endptr, 10);

	if(*endptr != 0){
		/* JP-182 */
		return 0;
	}

	return 1;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-183 */


static
int iff_field_to_int( char* field, int default_val, int* val )
{
	char* endptr;

	if( field == NULL || (field && *field == 0)) {
		*val = default_val;
		return 2;
	}

	*val = strtol( field, &endptr, 10);

	if(*endptr != 0){
		/* JP-184 */
		return 0;
	}

	return 1;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-185 */

static
int iff_field_to_double( char* field, double default_val, double* val )
{
	char* endptr;
	char buf[IFF_FIELD_SIZE];
	char *p, *bp;

	if( field == NULL || (field && *field == 0)) {
		*val = default_val;
		return 2;
	}

	*val = strtod( field, &endptr);

	if(*endptr != 0){
		/* JP-186 */
		
		p = field;
		bp = buf;
		while(*p){
			if( p == endptr ){
				if(*p == '+' || *p == '-'){
					*bp = 'E';
					bp++;
				} else
					return 0;
			}
			*bp = *p;
			p++;
			bp++;
		}
		*bp = 0;
		*val = strtod( buf, &endptr);

		if(*endptr != 0){
			/* JP-187 */
			return 0;
		}
	}

	return 1;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-188 */

static
int iff_field_to_param( char* field, char type, void* param, void** next_param_pos )
{
	int R;
	int* p_int;
	double* p_double;

	type = toupper(type);

	switch( type ){
	case 'I':
		R = iff_field_to_int( field, *((int*)param), (int*)param );
		p_int = (int*)param;
		p_int++;
		param = p_int;
		break;
	case 'U':
		R = iff_field_to_uint( field, *((unsigned int*)param), (unsigned int*)param );
		p_int = (int*)param;
		p_int++;
		param = p_int;
		break;
	case 'D':
		R = iff_field_to_double( field, *((double*)param), (double*)param );
		p_double = (double*)param;
		p_double++;
		param = p_double;
		break;
	case 'S':
		if( field == NULL || field[0] == 0)
			R = 2;
		else {
			strcpy( (char*)param, field );
			R = 1;
		}
	case '_':
		/* JP-189 */
		R = 1; /* JP-190 */
		break;
	default:
		HECMW_assert(0);
	}

	if( next_param_pos ){
		*next_param_pos = param;
	}

	return R;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-191 */

static
void iff_dump( char* iff )
{
	int i,j;
	char* p;

	putc('|',stdout);

	p = iff;
	for( i=0; i<FIELD_NUMBER; i++){
		for( j=0; j<IFF_FIELD_SIZE; j++){
			if(*p == 0)
				putc('_',stdout);
			else
				putc( *p, stdout);
			p++;
		}
		putc('|',stdout);
		p++;
	}
	putc('\n',stdout);
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-192 */

static
void iff_dump_csv( char* iff )
{
	int i;
	char *f;
	for( i=0; i<FIELD_NUMBER; i++){
		f = iff_get_field( iff, i );
		printf("%s,", f);
	}
	printf("\n");
}



/* iff_conv.c */ 

/* JP-193 */



/* JP-194 */



/* JP-195 */

static void free_to_iff_format( char* line, char* iff)
{
	/* JP-196 */

	char token[IFF_FIELD_BUFFER_SIZE];
	char *line_p, *token_p;
	int i;

	iff_clear(iff);

	i = 0;
	token_p = token;
	line_p = line;

	while( *line_p ) {
		if(*line_p == ','){
			*token_p = 0;
			if( i>=FIELD_NUMBER ){
				printf("フィールドが多すぎる");
				break;
			}
			iff_set_field( iff, i, token );
			token_p = token;
			i++;
		}else{
			*token_p = *line_p;
			token_p++;
		}
		line_p++;
	}

	if( i>0 && i<FIELD_NUMBER){
		*token_p = 0;
		iff_set_field( iff, i, token );
	}

	iff_regulize( iff );
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-197 */

static void small_to_iff_format( char* line, char* iff)
{
	int i, j;
	char token[SMALL_FIELD_COLUMN+1];

	iff_clear(iff);

	i = 0;

	if( get_fixed_token( line, token, FIRST_FIELD_COLUMN )){
		iff_set_field( iff, i, token );
	}else{
		/* JP-198 */
		printf("フィールドが多すぎる");
	}

	i++;

	for(;i<FIELD_NUMBER; i++){
		if( get_fixed_token( NULL, token, SMALL_FIELD_COLUMN ))
			iff_set_field( iff, i, token );
		else
			break;
	}

	iff_regulize( iff );
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-199 */


static void large_to_iff_format( char* line1, char* line2, char* iff)
{
	char token[LARGE_FIELD_COLUMN+1];

	int i = 0;
	int j;
	char buff[LARGE_FIELD_COLUMN+1];

	iff_clear(iff);

	/* line1 -------------- */


	/* JP-200 */

	if( get_fixed_token( line1, token, FIRST_FIELD_COLUMN )){
		iff_set_field( iff, i, token );
	}else{
		/* JP-201 */
	}

	i++;

	for(;i<LARGE_FIELD_LINE_FIELD_NUMBER-1; i++){
		if( get_fixed_token( NULL, token, LARGE_FIELD_COLUMN ))
			iff_set_field( iff, i, token );
		else
			break;
	}

	/* JP-202 */

	/* line2 -------------- */


	/* JP-203 */
	get_fixed_token( line2, token, FIRST_FIELD_COLUMN );

	i = LARGE_FIELD_LINE_FIELD_NUMBER-1;

	for(;i<FIELD_NUMBER; i++, j++){
		if( get_fixed_token( NULL, token, LARGE_FIELD_COLUMN ))
			iff_set_field( iff, i, token );
		else
			break;
	}

	/* JP-204 */

	if(i == FIELD_NUMBER-1){
		if( get_fixed_token( NULL, token, LAST_FIELD_COLUMN ))
			iff_set_field( iff, i, token );
	}

	iff_regulize( iff );/* JP-205 */
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-206 */

static 
char* get_fixed_token( char* src, char* token, int field_size )
{
	int j;
	char* token_p;
	static char* src_p = NULL;

	if( src )
		src_p = src;

	if( src_p == NULL || *src_p == 0 || *src_p == '\n')
		return NULL;

	token_p = token;

	for(j=0; j<field_size && *src_p != 0; j++, src_p++, token_p++){
		*token_p = *src_p;
	}
	*token_p = 0;

	return token;	
}


/* iff_node.c */ 

/* JP-207 */


/* JP-208 */

static iff_node_t* iff_node_create(void)
{
	iff_node_t * node = HECMW_malloc( sizeof( iff_node_t ) );
	iff_node_init(node);
	return node;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-209 */

static void iff_node_free( iff_node_t* node )
{
	HECMW_free( node );
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-210 */

static void iff_node_init( iff_node_t* node )
{
	iff_init(node->iff);
	node->line_no = -1;
	node->next = NULL;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-211 */

static void iff_node_set( iff_node_t* node, char* iff, int line_no)
{
	iff_copy( node->iff, iff );
	node->line_no = line_no;
}


/* read_iff.c */ 

/* JP-212 */
/* JP-213 */



#define read_iff_error -1
#define read_iff_comment 0
#define read_iff_success 1
#define read_iff_eof 2

static
int read_iff( FILE* fp, char* iff, int* line_no )
{
	char line1[ IFF_SIZE ];
	char line2[ IFF_SIZE ];
	
	char* fg = fgets( line1, IFF_SIZE-1, fp);
	if(!fg)
		return read_iff_eof;

	(*line_no)++;
	remove_cr(line1);

	switch( field_format_judge( line1 )) {
	case fft_free_field_format:	/* JP-214 */
		free_to_iff_format( line1, iff );
		break;
	case fft_small_field_format:	/* JP-215 */
		small_to_iff_format( line1, iff );
		break;
	case fft_large_field_format:	/* JP-216 */
		/* JP-217 */
		if(!fgets( line2, IFF_SIZE-1, fp)){
			return read_iff_error;	/* JP-218 */
		}
		(*line_no)++;
		remove_cr(line2);
		large_to_iff_format( line1, line2, iff );
		break;
	case fft_comment:
		return read_iff_comment;	/* JP-219 */
	default:
		return read_iff_error; /* JP-220 */
	}

	return read_iff_success;
}

/* bulk.c */ 

/* JP-221 */


/* JP-222 */

static void iff_bulk_init( iff_bulk_t* bulk )
{
	bulk->iff_node_root = NULL;
	bulk->iff_node_last = NULL;
	bulk->next_bulk = NULL;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-223 */

static iff_bulk_t* iff_bulk_create(void)
{
	iff_bulk_t* bulk = (iff_bulk_t*)HECMW_malloc( sizeof( iff_bulk_t ));
	iff_bulk_init( bulk );

	return bulk;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-224 */

static void iff_bulk_free( iff_bulk_t* bulk )
{
	iff_node_t* next;
	iff_node_t* node = bulk->iff_node_root;

	while( node ) {
		next = node->next;
		iff_node_free( node );
		node = next;
	}

	HECMW_free( bulk );
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-225 */
/* JP-226 */

static void iff_bulk_regist( iff_bulk_t* bulk, char* iff, int line_no )
{
	iff_node_t *node, *last;
	iff_node_t* new_node;

	new_node = iff_node_create();
	iff_node_set( new_node, iff, line_no);

	/* JP-227 */

	if(bulk->iff_node_last){
		bulk->iff_node_last->next = new_node;
		bulk->iff_node_last = new_node;
	}else{
		bulk->iff_node_last = bulk->iff_node_root = new_node;
	}
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-228 */

static void iff_bulk_move_append( iff_bulk_t* dest, iff_bulk_t* src )
{
	if(src->iff_node_root == NULL)
		return;

	if(dest->iff_node_last)
		dest->iff_node_last->next = src->iff_node_root;
	else
		dest->iff_node_root = src->iff_node_root;
	
	dest->iff_node_last = src->iff_node_last;

	src->iff_node_root = NULL;
	src->iff_node_last = NULL;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-229 */

static int iff_bulk_is_completed( iff_bulk_t* bulk )
{
	int fg;

	if( bulk->iff_node_root == NULL )
		return 0; /* JP-230 */

	fg = iff_is_parent( bulk->iff_node_root->iff )
		&& iff_is_last( bulk->iff_node_last->iff);

	return fg;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-231 */

static char* iff_bulk_get_name( iff_bulk_t* bulk )
{
	if(!bulk->iff_node_root)
		return NULL;

	return iff_get_field( bulk->iff_node_root->iff, 0);
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-232 */

static char* iff_bulk_get_field( iff_bulk_t* bulk, int field_no )
{
	int i, line, n;
	iff_node_t* node;

	line = field_no / FIELD_NUMBER;
	n = field_no % FIELD_NUMBER;

	i = 0;
	node = bulk->iff_node_root;
	while( node && i<line ) {
		node = node->next;
		i++;
	}
	if(i<line)
		return NULL;

	return iff_get_field( node->iff, n);
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-233 */
/* JP-234 */

static char* iff_bulk_get_param( iff_bulk_t* bulk, int param_no )
{
	int i;
	iff_node_t* node;
	int N = FIELD_NUMBER - 2;
	int line = param_no / N;
	int n = param_no % N + 1;

	i = 0;
	node = bulk->iff_node_root;
	while( node && i<line ) {
		node = node->next;
		i++;
	}
	if(node == NULL)
		return NULL;

	return iff_get_field( node->iff, n);
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-235 */
/* JP-236 */


static char* iff_bulk_get_param_data( iff_bulk_t* bulk, int param_no, int* line_no, int* field_no )
{
	int i;
	iff_node_t* node;
	int N = FIELD_NUMBER - 2;
	int line = param_no / N;
	int n = param_no % N + 1;

	i = 0;
	node = bulk->iff_node_root;
	while( node && i<line ) {
		node = node->next;
		i++;
	}
	if(node == NULL)
		return 0;

	*line_no = node->line_no;

	*field_no = n + 1;

	return iff_get_field( node->iff, n);
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-237 */
/* JP-238 */



/* JP-239 */

static int get_number_in_str( char** pp)
{
	char* p = *pp;
	char s[40];
	char* sp = s;
	int n;

	while( *p >= '0' && *p <= '9'){
		*sp = *p;
		sp++;
		p++;
	}
	*sp = 0;
	if( s[0] != 0)
		n = atoi(s);
	else
		n = 0;

	*pp = p;
	return n;
}


static int iff_bulk_get_param_list( iff_bulk_t* bulk, const char* format, int* result, ... )
{
	va_list ap;

	char* nf;
	int i;
	int arr_n;	/* JP-240 */
	int ness_n;	/* JP-241 */
	int arr_fg;	/* JP-242 */
	int R;
	int* r;
	char* f;
	char* field;
	char type;
	int ness;	/* JP-243 */

	int line_no;	/* JP-244 */
	int field_no;	/* JP-245 */

	int param_pos;
	void* param;
	char* param_name_p;
	char param_name[IFF_FIELD_SIZE];
	int read_param_count;

	va_start(ap, result);

	f = (char*)format;
	r = result;
	param_pos = 0;
	read_param_count = 0;

	while(*f){
		type = toupper(*f);
		if( type != '_' )
			ness = isupper(*f);
		else
			ness = 0;

		/* JP-246 */

		if( ness )
			ness_n = 1;
		else
			ness_n = 0;

		nf = f;
		nf++;
		arr_n = get_number_in_str( &nf );
		arr_fg = 1;
		if(arr_n==0){
			arr_n = 1;
			arr_fg = 0;
		} else if( *nf == '/'){
			nf++;
			ness_n = get_number_in_str( &nf );
			if(!ness)
				ness_n = 0;
		} else if( ness ){
			ness_n = arr_n;
		}

		HECMW_assert( !((type == 'S') && (arr_n > 1)) );

		/* JP-247 */

		if(type != '_') {
			param = va_arg(ap, void* );
			param_name_p = va_arg(ap, char* );
		}

		for(i=0; i<arr_n; i++){
			field = iff_bulk_get_param_data( bulk, param_pos, &line_no, &field_no );
			param_pos++;
			R = iff_field_to_param( field, type, param, &param );

			if( R == 0 ){				/* JP-248 */
				if( arr_fg )
					sprintf( param_name, "%s%d", param_name_p, i+1 );
				else
					strcpy( param_name, param_name_p);
				set_error_of_field_convert( grid_filename, line_no, field_no,
					iff_bulk_get_name(bulk), param_name, field, type );
				return -1;
			} else if( i<ness_n && R == 2 ){	/* JP-249 */
				if( arr_fg )
					sprintf( param_name, "%s%d", param_name_p, i+1 );
				else
					strcpy( param_name, param_name_p);
				set_error_of_blank_field( grid_filename, 
					line_no,  field_no, iff_bulk_get_name(bulk), param_name );
				return -1;
			}

			if( type != '_' ){
				*r = R;
				r++;
			}
		}

		f = nf;
	}

	va_end(ap);
	return 0;
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-250 */
/* JP-251 */





static int iff_bulk_get_param_list_pattern( iff_bulk_t* bulk, const char* format, int start_param_pos, int* size, ... )
{
	#define GPLP_PARAM_LIST_SIZE 50

	va_list ap;

	int i;
	char* field;
	char type;

	int param_pos;
	void* param;
	char param_name[IFF_FIELD_SIZE];

	void* param_list_p[ GPLP_PARAM_LIST_SIZE ];
	void** param_list[ GPLP_PARAM_LIST_SIZE ];
	char* param_name_head[ GPLP_PARAM_LIST_SIZE ];

	int R[ GPLP_PARAM_LIST_SIZE ];

	int param_list_size;	/* JP-252 */
	int data_number;	/* JP-253 */
	int data_size;		/* JP-254 */

	int index;
	int loop_fg;

	int type_size;
	int r;
	int line_no;
	int field_no;
	int state; /* JP-255 */

	/* ------------------------------------------- */

	va_start(ap, size);

	str_upr( (char*)format );
	param_list_size = strlen((char*)format);

	/* JP-256 */

	for(i=0; i<param_list_size; i++){
		param_list[i] = va_arg(ap, void** );
		param_name_head[i] = va_arg(ap, char* );
	}

	/* JP-257 */

	data_size = iff_bulk_get_line_number( bulk ) / param_list_size;
	for(i=0; i<param_list_size; i++){
		switch( format[i] ){
		case 'I':
			type_size = sizeof(int);
			break;
		case 'U':
			type_size = sizeof( unsigned int);
			break;
		case 'D':
			type_size = sizeof(double);
			break;
		default:
			HECMW_assert(0);
		}
		*param_list[i] = HECMW_malloc( type_size * data_size );
		param_list_p[i] = *param_list[i];
	}

	param_pos = start_param_pos;
	index = 0;
	loop_fg = 1;
	while(loop_fg){
		state = 0;
		for(i=0; i<param_list_size; i++){
			field = iff_bulk_get_param_data( bulk, param_pos, &line_no, &field_no );
			if( field == NULL || field[0] == 0){
				goto END;
				break;
			}
			param_pos++;
			type = format[i];
			param = param_list_p[i];
			r = iff_field_to_param( field, type, param, &param );
			if( r == 1 ){
				param_list_p[i] = param; /* JP-258 */
			} else if( r == 0){		/* error */
				sprintf( param_name, "%s%d", param_name_head[i], index+1 );
				set_error_of_field_convert( grid_filename, line_no, field_no,
							iff_bulk_get_name(bulk), param_name, field, type );
				return -1;
			}
			if( state != 0 && state != r ){
				sprintf( param_name, "%s%d", param_name_head[i], index+1 );
				set_error_of_field_convert( grid_filename, line_no, field_no,
							iff_bulk_get_name(bulk), param_name, field, type );
				return -1;
			}
			state = r;
		}
		index++;
	}
	END:
	*size = index;

	va_end(ap);
	return 0;

	#undef GPLP_PARAM_LIST_SIZE
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-259 */

static int iff_bulk_get_line_number( iff_bulk_t* bulk )
{
	int i = 0;
	iff_node_t* node = bulk->iff_node_root;

	while( node ) {
		i++;
		node = node->next;
	}

	return i;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-260 */

static int iff_bulk_get_line_no( iff_bulk_t* bulk, int iff_pos )
{
	int i = 0;
	iff_node_t* node = bulk->iff_node_root;

	while( node ) {
		if( i == iff_pos )
			return node->line_no;
		i++;
		node = node->next;
	}
	return 0;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-261 */

static iff_node_t* iff_bulk_search_pointing_line( iff_bulk_t* bulk, char* line_pointer )
{
	/* JP-262 */

	char *p;

	iff_node_t* node;

	node = bulk->iff_node_last;

	if( node && iff_is_pointing_line( node->iff )){
		p = iff_get_pointing_line_pointer( node->iff );
		if(strcmp(p, line_pointer) == 0){
			return node;
		}
	}

	return NULL;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-263 */

static iff_node_t* iff_bulk_search_continuous_line( iff_bulk_t* bulk, char* line_pointer )
{
	/* JP-264 */

	char *p;
	
	iff_node_t* node = bulk->iff_node_root;

	if( node && iff_is_continuous_line( node->iff)){
		p = iff_get_continuous_line_pointer( node->iff );
		if( strcmp( p, line_pointer) == 0){
			return node;
		}
	}

	return NULL;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-265 */


static void iff_bulk_dump( iff_bulk_t* bulk)
{
	iff_node_t* node = bulk->iff_node_root;

	while( node ) {
		printf( "%3d", node->line_no);
		iff_dump( node->iff );
		node = node->next;
	}

}

/* bulk_parse_head.c */ 

/* JP-266 */

/* JP-267 */

/* ------------------------------------------------------------------------------------------------ */


#include <stdarg.h>
#include <stdlib.h>


#define GENERATE_CODE( name ) \
	static int iff_bulk_parse_##name( iff_bulk_t* bulk );

	GENERATE_CODE( GRID )
	GENERATE_CODE( MAT1 )
	GENERATE_CODE( MAT4 )
	GENERATE_CODE( PROD )
	GENERATE_CODE( PSHELL )
	GENERATE_CODE( PSOLID )
	GENERATE_CODE( CROD )
	GENERATE_CODE( CBAR )
	GENERATE_CODE( CBEAM )
	GENERATE_CODE( CTRIA3 )
	GENERATE_CODE( CTRIA6 )
	GENERATE_CODE( CQUAD4 )
	GENERATE_CODE( CQUAD8 )
	GENERATE_CODE( CTETRA )
	GENERATE_CODE( CPENTA )
	GENERATE_CODE( CHEXA )
	GENERATE_CODE( CTRIAX6 )
	GENERATE_CODE( INCLUDE )

#undef GENERATE_CODE


#define GENERATE_CODE( name ) \
	BULK_TYPE_##name,

enum {
	GENERATE_CODE( GRID )
	GENERATE_CODE( MAT1 )
	GENERATE_CODE( MAT4 )
	GENERATE_CODE( PROD )
	GENERATE_CODE( PSHELL )
	GENERATE_CODE( PSOLID )
	GENERATE_CODE( CROD )
	GENERATE_CODE( CBAR )
	GENERATE_CODE( CBEAM )
	GENERATE_CODE( CTRIA3 )
	GENERATE_CODE( CTRIA6 )
	GENERATE_CODE( CQUAD4 )
	GENERATE_CODE( CQUAD8 )
	GENERATE_CODE( CTETRA )
	GENERATE_CODE( CPENTA )
	GENERATE_CODE( CHEXA )
	GENERATE_CODE( CTRIAX6 )
	GENERATE_CODE( INCLUDE )
	BULK_TYPE_UNKNOWN
};
#undef GENERATE_CODE


#define BULK_TYPE_NUMBER   BULK_TYPE_UNKNOWN

#define GENERATE_CODE( name ) \
	#name ,

static char *bulk_type_name[] = {
	GENERATE_CODE( GRID )
	GENERATE_CODE( MAT1 )
	GENERATE_CODE( MAT4 )
	GENERATE_CODE( PROD )
	GENERATE_CODE( PSHELL )
	GENERATE_CODE( PSOLID )
	GENERATE_CODE( CROD )
	GENERATE_CODE( CBAR )
	GENERATE_CODE( CBEAM )
	GENERATE_CODE( CTRIA3 )
	GENERATE_CODE( CTRIA6 )
	GENERATE_CODE( CQUAD4 )
	GENERATE_CODE( CQUAD8 )
	GENERATE_CODE( CTETRA )
	GENERATE_CODE( CPENTA )
	GENERATE_CODE( CHEXA )
	GENERATE_CODE( CTRIAX6 )
	GENERATE_CODE( INCLUDE )
	"UNKNOWN"
};
#undef GENERATE_CODE




/*-------------------------------------------------------*/


static int iff_bulk_type( const char* bulkname )
{
	int i;

	for(i=0; i<BULK_TYPE_NUMBER; i++){
		if( 0 == strcmp( bulkname, bulk_type_name[i] ))
			return i;
	}

	return 	BULK_TYPE_UNKNOWN;
}




static
int iff_bulk_parse( iff_bulk_t* bulk )
{
	char* bulkname;
	int no;

	bulkname = iff_bulk_get_name( bulk );

	HECMW_assert(bulkname);

	no = iff_bulk_type( bulkname );

	#define GENERATE_CODE( name ) \
		case BULK_TYPE_##name:\
			return iff_bulk_parse_##name( bulk );

	switch( no ){
		GENERATE_CODE( GRID )
		GENERATE_CODE( MAT1 )
		GENERATE_CODE( MAT4 )
		GENERATE_CODE( PROD )
		GENERATE_CODE( PSHELL )
		GENERATE_CODE( PSOLID )
		GENERATE_CODE( CROD )
		GENERATE_CODE( CBAR )
		GENERATE_CODE( CBEAM )
		GENERATE_CODE( CTRIA3 )
		GENERATE_CODE( CTRIA6 )
		GENERATE_CODE( CQUAD4 )
		GENERATE_CODE( CQUAD8 )
		GENERATE_CODE( CTETRA )
		GENERATE_CODE( CPENTA )
		GENERATE_CODE( CHEXA )
		GENERATE_CODE( CTRIAX6 )
		GENERATE_CODE( INCLUDE )
	default:
		/* JP-268 */
		/* JP-269 */
		return 0;
		/* return iff_bulk_parse_debug( bulk ); */
	}
}

/* bulk_parse.c */ 


/*---------------------------------------------------------------*/
/* JP-270 */

typedef struct ST_PSHELL_MID2 {
	int PID;
	int MID2;
	struct ST_PSHELL_MID2* next;
} pshell_mid2_t;


typedef struct ST_PSHELLMID2_LIST {
	pshell_mid2_t* root;
} pshell_mid2_list_t;


static pshell_mid2_list_t pshell_mid2_list;
static int pshell_mid2_list_init_fg = 0;

static void pshell_mid2_list_init(void);
static void pshell_mid2_list_clear(void);
static void pshell_mid2_list_regist( int pid, int mid2 );
static int pshell_mid2_list_get_mid2_by_pid( int pid, int* mid2 );


/* JP-271 */

static void pshell_mid2_list_init(void)
{
	pshell_mid2_list.root = NULL;
	pshell_mid2_list_init_fg = 1;
}

static void pshell_mid2_list_clear(void)
{
	pshell_mid2_t* node = pshell_mid2_list.root;
	pshell_mid2_t* next;

	while( node ){
		next = node->next;
		HECMW_free( node );
		node = next;
	}
	pshell_mid2_list.root = NULL;
}


static void pshell_mid2_list_regist( int pid, int mid2 )
{
	pshell_mid2_t* node;

	node = HECMW_malloc( sizeof( pshell_mid2_t ));
	HECMW_assert( node );

	node->PID = pid;
	node->MID2 = mid2;

	if(!pshell_mid2_list_init_fg){
		pshell_mid2_list_init();
	}


	if(pshell_mid2_list.root){
		node->next = pshell_mid2_list.root;
	}else{
		node->next = NULL;
	}
	pshell_mid2_list.root = node;
}


static int pshell_mid2_list_get_mid2_by_pid( int pid, int* mid2 )
{
	pshell_mid2_t* node = pshell_mid2_list.root;

	while( node ){
		if(node->PID == pid){
			*mid2 = node->MID2;
			return 1;
		}
		node = node->next;
	}

	return 0;
}


/* JP-272 */
/* JP-273 */

	
static int check_order_of_elem( int result[], int G_start_pos, int G_number, int G_ness )
{
	int i;
	int first_order = 0;
	int *r = &(result[ G_start_pos ]);

	for(i = 0; i<G_ness; i++, r++){
		if( *r != 1 ){
			return 0;	/* JP-274 */
		}
	}

	for(; i<G_number; i++, r++){
		if( *r != 1 ){
			return 1;	/* JP-275 */
		}
	}

	return 2;
}



/* JP-276 */
/* JP-277 */

/* ------------------------------------------------------------------------------------------------ */
/* JP-278 */



static
struct hecmw_io_material * create_mat_struct(const char* name, int item_n, int sub_n[], double data[] )
{
	#define cms_malloc( p, size ) { \
		(p) = HECMW_malloc((size));  \
		if((p) == NULL) { \
			HECMW_set_error(errno, ""); return NULL; \
		}\
	}

	struct hecmw_io_material *mat;
	struct hecmw_io_matitem *matitem;
	struct hecmw_io_matsubitem *matsubitem;
	int i,j, d_count;
	int val_n;

	cms_malloc( mat, sizeof(*mat));
	cms_malloc( matitem, sizeof(*matitem) * item_n );

	d_count = 0;

	for(i=0; i<item_n; i++){
		matitem[i].item = i+1;
		matitem[i].nval = sub_n[i];

		cms_malloc( matsubitem, sizeof(*matsubitem) ); 
		cms_malloc( matsubitem->val, sizeof(*(matsubitem->val)) * sub_n[i] );

		for(j=0; j<sub_n[i]; j++){
			matsubitem->val[j] = data[d_count];
			d_count++;
		}
		matsubitem->temp = data[d_count];
		d_count++;

		matsubitem->next = NULL;
		matitem[i].subitem = matsubitem;
	}

	strcpy( mat->name, name);
	mat->nitem = item_n;
	mat->item = matitem;
	mat->next = NULL;

	return mat;

	#undef cms_malloc
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-279 */


static
int iff_bulk_parse_MAT1( iff_bulk_t* bulk )
{
	#define BULK_PARAM_NUMBER 12

	unsigned int MID;	/* JP-280 */
	double E;		/* JP-281 */
	double G;		/* JP-282 */
	double NU;		/* JP-283 */
	double RHO;		/* JP-284 */
	double A;		/* JP-285 */
	double TREF;		/* JP-286 */
	double GE;		/* JP-287 */
	double ST,SC,SS;	/* JP-288 */
	unsigned int MCSID;	/* JP-289 */

	char format[] = "UDdDDDdddddu";

	int result[ BULK_PARAM_NUMBER ];

	MID = 0;
	E = G = NU = RHO = A = TREF = GE = ST = SC = SS = 0.0;
	MCSID = 0;

	iff_bulk_get_param_list( bulk, format, result, &MID, "MID", &E, "E", &G, "G", &NU, "NU",
							&RHO, "RHO", &A, "A", &TREF, "TREF", &GE, "GE",
							&ST, "ST", &SC, "SC", &SS,"SS", &MCSID, "MCSID");

	{
		int i;
		struct hecmw_io_material *mat;
		char name[HECMW_NAME_LEN+1];
		#define item_n 3
		int sub_n[item_n] = {2,1,1};
		double data[4 + item_n];

		matrial_name_by_MID( MID, name);

		i = 0;
		data[i++] = E;
		data[i++] = NU;
		data[i++] = 0; /* temp */
		data[i++] = RHO;
		data[i++] = 0; /* temp */
		data[i++] = A;
		data[i++] = 0; /* temp */

		mat = create_mat_struct( name, item_n, sub_n, data );
		if(HECMW_io_add_mat(mat) == NULL) return -1;

		#undef item_n
	}

	return 0;
	#undef BULK_PARAM_NUMBER
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-290 */


static
int iff_bulk_parse_MAT4( iff_bulk_t* bulk )
{
	#define BULK_PARAM_NUMBER 11

	unsigned int MID;	/* JP-291 */
	double K=0;		/* JP-292 */
	double CP=0;		/* JP-293 */
	double r=1.0;		/* JP-294 */
	double H;		/* JP-295 */
	double m;		/* JP-296 */
	double HGEN;		/* JP-297 */
	double REFENTH;		/* JP-298 */
	double TCH;		/* JP-299 */
	double TDELTA;		/* JP-300 */
	double QLAT;		/* JP-301 */

	char format[] = "Udddddddddd";

	int result[ BULK_PARAM_NUMBER ];

	iff_bulk_get_param_list( bulk, format, result, &MID, "MID", &K, "K", &CP, "CP", &r, "rho",
		&H, "H", &m, "m", &HGEN, "HGEN", &REFENTH, "REFENTH", &TCH, "TCH", &TDELTA, "TDELTA", &QLAT, "QLAT");

	{
		int i;
		struct hecmw_io_material *mat;
		char name[HECMW_NAME_LEN+1];
		#define item_n 3
		int sub_n[item_n] = {1,1,1};
		double data[3 + item_n];

		matrial_name_by_MID( MID, name);

		i = 0;
		data[i++] = K;
		data[i++] = 0; /* temp */
		data[i++] = CP;
		data[i++] = 0; /* temp */
		data[i++] = r;
		data[i++] = 0; /* temp */

		mat = create_mat_struct( name, item_n, sub_n, data );
		if(HECMW_io_add_mat(mat) == NULL) return -1;

		#undef item_n
	}

	return 0;
	#undef BULK_PARAM_NUMBER
}





/* JP-302 */
/* JP-303 */

/* ------------------------------------------------------------------------------------------------ */
/* JP-304 */

/* JP-305 */

typedef struct hecmw_io_section section_rec ;


static
void init_section_rec( section_rec* a)
{
	memset( a, 0, sizeof( section_rec ));
}

CREATE_LIST_CODE( section_rec )

static
int is_eq_section_rec( section_rec a, section_rec b)
{
	return (strcmp(a.egrp, b.egrp) == 0) ;
}

static
void copy_section_rec( section_rec* dest, section_rec src)
{
	memcpy( dest, &src, sizeof( section_rec ));
}

static
void delete_section_rec( section_rec a )
{
	/* JP-306 */
}

static t_list_section_rec section_list;

/*---------------------------------------------*/

/* JP-307 */

static void section_list_init(void)
{
	list_init_section_rec( &section_list );
}


/* JP-308 */

static int section_list_regist( struct hecmw_io_section *sect )
{
	return ( list_add_section_rec( &section_list, *sect ) == NULL);
}


/* JP-309 */

static int section_list_finalize(void)
{
	list_scan_start_section_rec( &section_list );
	while( list_scan_next_section_rec( &section_list )){
		section_rec* rec = &section_list.current->data;
		struct hecmw_io_egrp *egrp = HECMW_io_get_egrp( rec->egrp );
		if( egrp ){
			if(HECMW_io_add_sect(rec) == NULL)
				return -1;
		}
	}

	list_clear_section_rec( &section_list );

	return 0;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-310 */


static int iff_bulk_parse_PROD( iff_bulk_t* bulk )
{
	#define BULK_PARAM_NUMBER 6

	unsigned int PID; 	/* JP-311 */
	unsigned int MID;	/* JP-312 */
	double A;		/* JP-313 */
	double J;		/* JP-314 */
	double C = 0;		/* JP-315 */
	double NSM;		/* JP-316 */

	char format[] = "UUdddd";
	int result[ BULK_PARAM_NUMBER ];
	char grp_name[HECMW_NAME_LEN+1];
	struct hecmw_io_section sect;

	iff_bulk_get_param_list( bulk, format, result,
				&PID, "PID", &MID, "MID", &A, "A", &J, "J", &C, "C", &NSM, "NSM" );

	strcpy( sect.egrp, egrp_name_by_PID( PID, grp_name ));
	strcpy( sect.material, matrial_name_by_MID( MID, grp_name ));
	sect.composite = 1;
	sect.secopt = 0;

	sect.type = HECMW_SECT_TYPE_SOLID;
	sect.sect.solid.thickness = A;

	sect.next = NULL;

	if(section_list_regist(&sect)) return -1;

	return 0;

	#undef BULK_PARAM_NUMBER
}





/* ------------------------------------------------------------------------------------------------ */
/* JP-317 */


static int iff_bulk_parse_PSHELL( iff_bulk_t* bulk )
{
	#define BULK_PARAM_NUMBER 11

	unsigned int PID; 	/* JP-318 */
	unsigned int MID1;	/* JP-319 */
	double T;	 	/* JP-320 */
	int MID2;		/* JP-321 */
	double _12I_T_3;	/* JP-322 */
	unsigned int MID3;	/* JP-323 */
	double TS_S;		/* JP-324 */
	double NSM;		/* JP-325 */
	double Z1,Z2;		/* JP-326 */
	unsigned int MID4;	/* JP-327 */

	char format[] = "Uudiduddddi";

	int result[ BULK_PARAM_NUMBER ];

	char grp_name[HECMW_NAME_LEN+1];

	struct hecmw_io_section sect;

	/* JP-328 */

	MID1 = 0;
	MID2 = 0;
	_12I_T_3 = 1.0;
	TS_S = 0.833333;

	iff_bulk_get_param_list( bulk, format, result,
		&PID, "PID", &MID1, "MID1", &T, "T", &MID2, "MID2", &_12I_T_3, "12I/T^3", &MID3, "MID3", &TS_S,"TS/S", &NSM, "NSM",
		&Z1, "Z1", &Z2, "Z2", &MID4, "MID4" );

	strcpy( sect.egrp, egrp_name_by_PID( PID, grp_name ));
	strcpy( sect.material, matrial_name_by_MID( MID1, grp_name ));
	sect.composite = 1;
	sect.secopt = 0;

	if( MID2 == 0 )
		sect.secopt = HECMW_SECT_OPT_PSTRESS;
	else if( MID2 == -1 )
		sect.secopt = HECMW_SECT_OPT_PSTRAIN;

	if( MID2 == 0 || MID2 == -1)
		sect.type = HECMW_SECT_TYPE_SOLID;
	else
		sect.type = HECMW_SECT_TYPE_SHELL;

	if(sect.type == HECMW_SECT_TYPE_SOLID){
		sect.sect.solid.thickness = T;
	} else {
		sect.sect.shell.thickness = T;
		sect.sect.shell.integpoints = 1; /* JP-329 */
	}

	sect.next = NULL;

	if(section_list_regist(&sect)) return -1;

	/* JP-330 */

	pshell_mid2_list_regist( PID,  MID2);

	return 0;

	#undef BULK_PARAM_NUMBER
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-331 */


static int iff_bulk_parse_PSOLID( iff_bulk_t* bulk )
{
	#define BULK_PARAM_NUMBER 7

	unsigned int PID;		/* JP-332 */
	unsigned int MID;		/* JP-333 */
	int CORDM;			/* JP-334 */
	char IN[IFF_FIELD_SIZE];	/* JP-335 */
	char STRESS[IFF_FIELD_SIZE];	/* JP-336 */
	char ISOP[IFF_FIELD_SIZE];	/* JP-337 */
	char FCTN[IFF_FIELD_SIZE];	/* JP-338 */


	char format[] = "UUissss";
	int result[ BULK_PARAM_NUMBER ];
	struct hecmw_io_section sect;
	char grp_name[HECMW_NAME_LEN+1];

	/* JP-339 */

	CORDM = 0;
	IN[0] = 0;
	STRESS[0] = 0;
	ISOP[0] = 0;
	strcpy( FCTN, "SMECH" );

	iff_bulk_get_param_list( bulk, format, result,
		&PID, "PID", &MID, "MID", &CORDM, "CORDM", IN, "IN", STRESS, "STRESS", ISOP, "ISOP", FCTN, "FCTN" );

	sect.type = HECMW_SECT_TYPE_SOLID;

	strcpy( sect.egrp, egrp_name_by_PID( PID, grp_name ));
	strcpy( sect.material, matrial_name_by_MID( MID, grp_name ));
	sect.composite = 1; /* JP-340 */

	sect.secopt = 0;  /* JP-341 */

	sect.sect.solid.thickness = 1; /* JP-342 */

	sect.next = NULL;

	if(section_list_regist(&sect)) return -1;

	return 0;

	#undef BULK_PARAM_NUMBER
}





/* JP-343 */
/* JP-344 */

/* ------------------------------------------------------------------------------------------------ */
/* JP-345 */



static int iff_bulk_parse_CROD( iff_bulk_t* bulk )
{
	#define G_NUMBER 2
	#define BULK_PARAM_NUMBER (2 + G_NUMBER)

	unsigned int EID;	/* JP-346 */
	unsigned int PID;	/* JP-347 */
	unsigned int G[ G_NUMBER ];

	char format[20] = "UUU2";
	int result[ BULK_PARAM_NUMBER ];
	char grp_name[HECMW_NAME_LEN+1];
	int EID_List[1];

	iff_bulk_get_param_list( bulk, format, result, &EID, "EID", G, "G");

	if( HECMW_io_add_elem( EID, HECMW_ETYPE_ROD1, (int*)G, 0, NULL) == NULL)
		return -1;

	EID_List[0] = EID;
	if(HECMW_io_add_egrp( egrp_name_by_PID( PID, grp_name ), 1, EID_List) == -1)
		return -1;

	if(HECMW_io_add_egrp( "ALL", 1, EID_List) == -1)
		return -1;

	return 0;

	#undef G_NUMBER
	#undef BULK_PARAM_NUMBER
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-348 */


static int iff_bulk_parse_CBAR( iff_bulk_t* bulk )
{
	return -1;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-349 */



static int iff_bulk_parse_CBEAM( iff_bulk_t* bulk )
{
	return -1;
}



/* JP-350 */
/* JP-351 */

/* ------------------------------------------------------------------------------------------------ */
/* JP-352 */



/* JP-353 */



static
int surface_elem_store( int bulk_type, unsigned int EID, unsigned int PID,
				unsigned int G[], int G_number, double THETA_MCID, double ZOFFS, double T[], int T_number)
{
	int elem_type;
	struct hecmw_io_egrp *egrp;
	char grp_name[HECMW_NAME_LEN+1];
	int EID_List[1];

	/* JP-354 */

	switch( bulk_type ){
	case BULK_TYPE_CTRIA3:
		elem_type = HECMW_ETYPE_TRI1;
		break;
	case BULK_TYPE_CTRIA6:
		elem_type = HECMW_ETYPE_TRI2;
		break;
	case BULK_TYPE_CQUAD4:
		elem_type = HECMW_ETYPE_QUA1;
		break;
	case BULK_TYPE_CQUAD8:
		elem_type = HECMW_ETYPE_QUA2;
		break;
	default:
		HECMW_assert(0);
	}

	if( HECMW_io_add_elem( EID, elem_type, (int*)G, T_number, T) == NULL)
		return -1;

	EID_List[0] = EID;

	if(HECMW_io_add_egrp( "ALL", 1, EID_List) == -1)
		return -1;

	if(HECMW_io_add_egrp( egrp_name_by_PID( PID, grp_name ), 1, EID_List) == -1)
		return -1;

	return 0;
}


/* JP-355 */
/* JP-356 */
/* JP-357 */




static
int surface_elem_type_decide(void)
{
	int i;
	struct hecmw_io_id_array *id_array;
	char grp_name[HECMW_NAME_LEN+1];
	pshell_mid2_t* node;
	struct hecmw_io_element* elem;

	node = pshell_mid2_list.root;
	while( node ){
		egrp_name_by_PID( node->PID, grp_name );
		id_array = HECMW_io_get_elem_in_egrp( grp_name );
		if(id_array){
			for(i=0; i<id_array->n; i++){
				elem = HECMW_io_get_elem( id_array->id[i] );
				if(!elem)
					HECMW_assert(0);
				if( node->MID2 >= 1 ){
					switch( elem->type ){
					case HECMW_ETYPE_TRI1:
							elem->type = HECMW_ETYPE_SHT1;
							break;
					case HECMW_ETYPE_TRI2:
							elem->type = HECMW_ETYPE_SHT2;
							break;
					case HECMW_ETYPE_QUA1:
							elem->type = HECMW_ETYPE_SHQ1;
							break;
					case HECMW_ETYPE_QUA2:
							elem->type = HECMW_ETYPE_SHQ2;
							break;
					}
				}
			}
			HECMW_free(id_array);
		}
		node = node->next;
	}

	pshell_mid2_list_clear();
	return 0;
}




/* ------------------------------------------------------------------------------------------------ */
/* JP-358 */


static int iff_bulk_parse_CTRIA3( iff_bulk_t* bulk )
{
	int i;

	#define BULK_PARAM_NUMBER (2+3+2+3)

	unsigned int EID;	/* JP-359 */
	unsigned int PID;	/* JP-360 */
	unsigned int G[3];	/* JP-361 */
	double THETA_MCID;	/* JP-362 */
	double ZOFFS;		/* JP-363 */
	double T[3];		/* JP-364 */
	int Tnumber;

	char format[] = "UUU3dd_3d3";
	int result[ BULK_PARAM_NUMBER ];

	T[0] = T[1] = T[2] = 0.0;

	iff_bulk_get_param_list( bulk, format, result, 
				&EID, "EID", &PID, "PID", G, "G", &THETA_MCID, "THETA/MCID", &ZOFFS, "ZOFFS", T, "T" );

	/* JP-365 */

	if( result[BULK_PARAM_NUMBER-3] == 1 && result[BULK_PARAM_NUMBER-2] == 1 && result[BULK_PARAM_NUMBER-1] == 1)
		Tnumber = 3;
	else
		Tnumber = 0;

	if(surface_elem_store( BULK_TYPE_CTRIA3, EID, PID, G, 3, THETA_MCID, ZOFFS, T, Tnumber))
		return -1;

	return 0;

	#undef BULK_PARAM_NUMBER
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-366 */



static int iff_bulk_parse_CTRIA6( iff_bulk_t* bulk )
{
	int i;

	#define G_NUMBER 6
	#define BULK_PARAM_NUMBER (2+G_NUMBER+2+3)

	unsigned int EID;		/* JP-367 */
	unsigned int PID;		/* JP-368 */
	unsigned int G[G_NUMBER];	/* JP-369 */
	double THETA_MCID;		/* JP-370 */
	double ZOFFS;			/* JP-371 */
	double T[3];			/* JP-372 */
	int Tnumber;
	int order;
	int conv_table[] = { 1, 2, 3, 5, 6, 4 };
	unsigned int node[G_NUMBER];

	char format[] = "UUU6dd_3d3";
	int result[ BULK_PARAM_NUMBER ];

	T[0] = T[1] = T[2] = 0.0;

	iff_bulk_get_param_list( bulk, format, result, 
				&EID, "EID", &PID, "PID", G, "G", &THETA_MCID, "THETA/MCID", &ZOFFS, "ZOFFS", T, "T" );

	/* JP-373 */

	if( result[BULK_PARAM_NUMBER-3] == 1 && result[BULK_PARAM_NUMBER-2] == 1 && result[BULK_PARAM_NUMBER-1] == 1)
		Tnumber = 3;
	else
		Tnumber = 0;


	for(i=0; i<G_NUMBER; i++)
		node[i] = G[ conv_table[i]-1];

	if(surface_elem_store( BULK_TYPE_CTRIA6, EID, PID, node, G_NUMBER, THETA_MCID, ZOFFS, T, Tnumber))
		return -1;

	return 0;

	#undef BULK_PARAM_NUMBER
	#undef G_NUMBER
}


/* ------------------------------------------------------------------------------------------------ */

/* JP-374 */



static int iff_bulk_parse_CQUAD4( iff_bulk_t* bulk )
{
	int i;

	#define G_NUMBER 4
	#define T_NUMBER 4
	#define BULK_PARAM_NUMBER (2+G_NUMBER+2+T_NUMBER)

	unsigned int EID;		/* JP-375 */
	unsigned int PID;		/* JP-376 */
	unsigned int G[G_NUMBER];	/* JP-377 */
	double THETA_MCID;		/* JP-378 */
	double ZOFFS;			/* JP-379 */
	double T[T_NUMBER];		/* JP-380 */
	int Tnumber;
	int order;

	char format[] = "UUU4dd_2d4";
	int result[ BULK_PARAM_NUMBER ];

	T[0] = T[1] = T[2] = T[3] = 0.0;

	iff_bulk_get_param_list( bulk, format, result, 
				&EID, "EID", &PID, "PID", G, "G", &THETA_MCID, "THETA/MCID", &ZOFFS, "ZOFFS", T, "T" );

	/* JP-381 */

	if( result[BULK_PARAM_NUMBER-4] == 1 && result[BULK_PARAM_NUMBER-3] == 1
	&& result[BULK_PARAM_NUMBER-2] == 1 && result[BULK_PARAM_NUMBER-1] == 1)
		Tnumber = T_NUMBER;
	else
		Tnumber = 0;

	if(surface_elem_store( BULK_TYPE_CQUAD4, EID, PID, G, G_NUMBER, THETA_MCID, ZOFFS, T, Tnumber))
		return -1;

	return 0;

	#undef BULK_PARAM_NUMBER
	#undef G_NUMBER
	#undef T_NUMBER
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-382 */




static int iff_bulk_parse_CQUAD8( iff_bulk_t* bulk )
{
	int i;

	#define G_NUMBER 8
	#define T_NUMBER 4
	#define BULK_PARAM_NUMBER (2+G_NUMBER+T_NUMBER+2)

	unsigned int EID;		/* JP-383 */
	unsigned int PID;		/* JP-384 */
	unsigned int G[G_NUMBER];	/* JP-385 */
	double T[T_NUMBER];		/* JP-386 */
	double THETA_MCID;		/* JP-387 */
	double ZOFFS;			/* JP-388 */
	int Tnumber;
	int order;

	char format[] = "UUU8d4dd";
	int result[ BULK_PARAM_NUMBER ];

	T[0] = T[1] = T[2] = T[3] = 0.0;

	iff_bulk_get_param_list( bulk, format, result, 
				&EID, "EID", &PID, "PID", G, "G", T, "T", &THETA_MCID, "THETA/MCID", &ZOFFS, "ZOFFS" );

	/* JP-389 */

	if( result[BULK_PARAM_NUMBER-4] == 1 && result[BULK_PARAM_NUMBER-3] == 1
	&& result[BULK_PARAM_NUMBER-2] == 1 && result[BULK_PARAM_NUMBER-1] == 1)
		Tnumber = T_NUMBER;
	else
		Tnumber = 0;

	if(surface_elem_store( BULK_TYPE_CQUAD8, EID, PID, G, G_NUMBER, THETA_MCID, ZOFFS, T, Tnumber))
		return -1;

	return 0;

	#undef BULK_PARAM_NUMBER
	#undef G_NUMBER
	#undef T_NUMBER
}




/* JP-390 */
/*
*/
/* ------------------------------------------------------------------------------------------------ */
/* JP-391 */


static
int iff_bulk_parse_solid_elem( iff_bulk_t* bulk,
				int first_etype, int second_etype,
				int conv_table[],
				int g_number, int g_ness_number )
{
	#define G_LIST_SIZE 100
	#define BULK_PARAM_LIST_SIZE (2+G_LIST_SIZE)

	unsigned int EID;	/* JP-392 */
	unsigned int PID;	/* JP-393 */

	unsigned int G[ G_LIST_SIZE ];
	unsigned int node[ G_LIST_SIZE ];

	int order; /* JP-394 */
	int i;

	char format[20];
	int result[ BULK_PARAM_LIST_SIZE ];

	char grp_name[HECMW_NAME_LEN+1];
	int EID_List[1];

	sprintf(format, "UUU%d/%d", g_number, g_ness_number );

	iff_bulk_get_param_list( bulk, format, result, &EID, "EID", &PID, "PID", G, "G"  );

	order = check_order_of_elem( result, 2, g_number, g_ness_number );

	/* JP-395 */

	for(i = 0; i<g_number; i++)
		node[i] = G[conv_table[i]-1];

	if( order == 1 ){
		if( HECMW_io_add_elem( EID, first_etype, (int*)node, 0, NULL) == NULL)
			return -1;
	}else if( order == 2 ){
		if( HECMW_io_add_elem( EID, second_etype, (int*)node, 0, NULL) == NULL)
			return -1;
	}else
		return -1;

	/* JP-396 */

	EID_List[0] = EID;
	if(HECMW_io_add_egrp( egrp_name_by_PID( PID, grp_name ), 1, EID_List) == -1)
		return -1;

	if(HECMW_io_add_egrp( "ALL", 1, EID_List) == -1)
		return -1;

	return 0;

	#undef G_LIST_SIZE
	#undef BULK_PARAM_LIST_SIZE
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-397 */



static int iff_bulk_parse_CTETRA( iff_bulk_t* bulk )
{
	int conv_table[] = {
				1, 2, 3, 4, 6, 7, 5, 8, 9, 10
	};
	return iff_bulk_parse_solid_elem( bulk,	HECMW_ETYPE_TET1, HECMW_ETYPE_TET2, conv_table, 10, 4);
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-398 */



static int iff_bulk_parse_CPENTA( iff_bulk_t* bulk )
{
	int conv_table[] = {
			1, 2, 3, 4, 5, 6, 8, 9, 7, 14, 15, 13, 10, 11, 12
	};
	return iff_bulk_parse_solid_elem( bulk,	HECMW_ETYPE_PRI1, HECMW_ETYPE_PRI2, conv_table, 15, 6);
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-399 */



static int iff_bulk_parse_CHEXA( iff_bulk_t* bulk )
{
	int conv_table[] = {
		1, 2, 3, 4, 5, 6, 7, 8,  9,10,11,12, 17,18,19,20, 13,14,15,16
	};
	return iff_bulk_parse_solid_elem( bulk,	HECMW_ETYPE_HEX1, HECMW_ETYPE_HEX2, conv_table, 20, 8);
}



/* JP-400 */
/*
*/
/* ------------------------------------------------------------------------------------------------ */
/* JP-401 */




static
int iff_bulk_parse_CTRIAX6( iff_bulk_t* bulk )
{
	#define G_NUMBER 6
	#define BULK_PARAM_NUMBER (2 + G_NUMBER + 1)

	int i;
	int conv_table[] = { 1, 3, 5, 4, 6, 2 };

	unsigned int EID;	/* JP-402 */
	unsigned int MID;	/* JP-403 */
	unsigned int G[ G_NUMBER ];
	double TH = 0.0;	/* JP-404 */

	unsigned int node[ G_NUMBER ];

	char format[20] = "UUU6d";
	int result[ BULK_PARAM_NUMBER ];

	iff_bulk_get_param_list( bulk, format, result, &EID, "EID", &MID, "MID", G, "G", &TH, "TH"  );

	/* JP-405 */

	for(i = 0; i< G_NUMBER; i++)
		node[i] = G[conv_table[i]-1];

	if( HECMW_io_add_elem( EID, HECMW_ETYPE_TRI2, (int*)node, 0, NULL) == NULL)
		return -1;

	/* JP-406 */

	{
		struct hecmw_io_egrp* egrp;
		char grp_name[HECMW_NAME_LEN+1];
		int EID_List[1];

		egrp_CTRIAX6_name_by_MID( MID, grp_name );
		egrp = HECMW_io_get_egrp(grp_name);

		if( !egrp ){ /* JP-407 */
			char mat_name[HECMW_NAME_LEN+1];
			struct hecmw_io_section sect;
			strcpy( sect.egrp, grp_name );
			strcpy( sect.material, matrial_name_by_MID( MID, mat_name ));
			sect.composite = 1; /* JP-408 */
			sect.secopt = HECMW_SECT_OPT_ASYMMETRY; /* JP-409 */
			sect.type = HECMW_SECT_TYPE_SOLID;
			sect.sect.solid.thickness = 1; /* JP-410 */
			sect.next = NULL;
			if(HECMW_io_add_sect(&sect) == NULL)
				return -1;
		}

		/* JP-411 */

		EID_List[0] = EID;
		if(HECMW_io_add_egrp( grp_name , 1, EID_List) == -1)
			return -1;

		if(HECMW_io_add_egrp( "ALL", 1, EID_List) == -1)
			return -1;
	}
	return 0;

	#undef G_NUMBER
	#undef BULK_PARAM_NUMBER
}




/* JP-412 */


static int iff_bulk_parse_GRID( iff_bulk_t* bulk )
{
	#define BULK_PARAM_NUMBER 8

	unsigned int ID=0;	/* JP-413 */
	unsigned int CP=0;	/* JP-414 */
	double X[3];		/* JP-415 */
	unsigned int CDI=0;	/* JP-416 */
	unsigned int PS=0;	/* JP-417 */
	unsigned int SEID;	/* JP-418 */

	char format[] = "UUD3uuu";
	int result[ BULK_PARAM_NUMBER ];

	int ID_List[1];

	X[0] = X[1] = X[2] = 0.0;

	iff_bulk_get_param_list( bulk, format, result, &ID, "ID", &CP, "CP", X, "X", &CDI,"CDI", &PS, "PS", &SEID, "SEID" );

	/* JP-419 */

	if( HECMW_io_add_node( ID, X[0], X[1], X[2]) == NULL )		/* JP-420 */
		return -1;

	ID_List[0] = ID;
	if(HECMW_io_add_ngrp( "ALL", 1, ID_List) == -1)
		return -1;

	#undef BULK_PARAM_NUMBER
	return 0;
}




/* JP-421 */


/* ------------------------------------------------------------------------------------------------ */
/*
	Nastran Quick Reference P.703
*/


static int iff_bulk_parse_INCLUDE( iff_bulk_t* bulk )
{
	int i;
	char buff[ IFF_SIZE ];
	char* iff;
	iff_node_t* node;
	char filename[HECMW_FILENAME_LEN+1];
	char *p, *bp;
	int fg;
	int errcode;

	/* JP-422 */

	node = bulk->iff_node_root;
	HECMW_assert(node);

	iff = node->iff;

	buff[0] = 0;
	for(i=0; i<FIELD_NUMBER;i++){
		strcat( buff, iff_get_field( iff, i ));
	}

	/* JP-423 */

	bp = buff;
	p = filename;
	fg = 0;
	while(*bp){
		if(fg){
			if( *bp == '\''){
				*p = 0;
				goto NEXT;
			}
			*p = *bp;
			p++;
		} else if( *bp=='\''){
			fg = 1;
		}
		bp++;
	}

	/* JP-424 */
	return -1;

	NEXT:


	/* JP-425 */

	printf("Include %s ==================================\n", filename);

	errcode = nastran_file_open( filename );
	return errcode;
}



/* bulk_list.c */ 

/* JP-426 */




/* ------------------------------------------------------------------------------------------------ */
/* JP-427 */

static void iff_bulk_list_init( iff_bulk_list_t* bulk_list)
{
	bulk_list->bulk_root = NULL;
	bulk_list->bulk_number = 0;
}
	
/* JP-428 */

static void iff_bulk_list_clear( iff_bulk_list_t* bulk_list)
{
	iff_bulk_t *node = bulk_list->bulk_root;
	iff_bulk_t *next;

	while( node ){
		next = node->next_bulk;
		iff_bulk_free( node );
		node = next;
	}
	bulk_list->bulk_root = NULL;
	bulk_list->bulk_number = 0;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-429 */

static int iff_bulk_list_is_empty( iff_bulk_list_t* list)
{
	return (list->bulk_root == NULL);
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-430 */

static void iff_bulk_list_regist( iff_bulk_list_t* bulk_list, iff_bulk_t* bulk)
{
	if( bulk_list->bulk_root == NULL) {
		bulk_list->bulk_root = bulk;
	} else {
		bulk->next_bulk = bulk_list->bulk_root;
		bulk_list->bulk_root = bulk;
	}	
	bulk_list->bulk_number++;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-431 */

static iff_bulk_t* iff_bulk_list_search_pointing_bulk( iff_bulk_list_t* bulk_list, char* line_pointer)
{
	iff_bulk_t* node = bulk_list->bulk_root;

	while(node){
		if( iff_bulk_search_pointing_line( node, line_pointer ))
			return node;
		node = node->next_bulk;
	}

	return NULL;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-432 */

static iff_bulk_t* iff_bulk_list_search_continuous_bulk( iff_bulk_list_t* bulk_list, char* line_pointer)
{
	iff_bulk_t* node = bulk_list->bulk_root;

	while(node){
		if( iff_bulk_search_continuous_line( node, line_pointer ))
			return node;
		node = node->next_bulk;
	}

	return NULL;
}

/* ------------------------------------------------------------------------------------------------ */
/* JP-433 */
/* JP-434 */

static int iff_bulk_list_remove( iff_bulk_list_t* bulk_list, iff_bulk_t* bulk )
{
	iff_bulk_t *prev, *node;
	node = bulk_list->bulk_root;
	prev = NULL;

	while(node){
		if( node == bulk ){
			if(prev)
				prev->next_bulk = node->next_bulk;
			else {
				bulk_list->bulk_root = node->next_bulk;
			}
			node->next_bulk = NULL;
			bulk_list->bulk_number--;
			return 1;
		}
		prev = node;
		node = node->next_bulk;
	}

	return 0;
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-435 */

#define TRYAL_NUMBER 20

static int iff_bulk_list_after_operation( iff_bulk_list_t* bulk_list )
{
	int bulk_number = bulk_list->bulk_number;
	iff_bulk_t* bulk;
	int counter = 0;

	while( !iff_bulk_list_is_empty( bulk_list ) ){
		bulk = bulk_list->bulk_root;
		while( bulk ){
			iff_bulk_t* next_bulk = bulk->next_bulk;
			if( iff_bulk_is_completed( bulk )){
				if( iff_bulk_parse( bulk )){
					iff_bulk_list_remove( bulk_list, bulk );
					iff_bulk_free( bulk );
				}
			} else {
				/* JP-436 */ /* JP-437 */
				return -1;
			}
			bulk = next_bulk;
		}
		if( bulk_number == bulk_list->bulk_number){ /* JP-438 */
			/* JP-439 */
			return -1;
		}
		bulk_number = bulk_list->bulk_number;
		counter++;
		if( counter == TRYAL_NUMBER ){
			/* JP-440 */
			return -1;
		}
	}

	return 0;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-441 */

static void iff_bulk_list_dump( iff_bulk_list_t * bulk_list )
{
	iff_bulk_t* node = bulk_list->bulk_root;

	printf(" Bulk List Dump ==============\n");

	while(node){
		iff_bulk_dump(node);
		node = node->next_bulk;
		printf("==============\n");
	}
}

/* iff_operation.c */ 

/* JP-442 */

/* JP-443 */

/* ================================================================================================================= */


static
int iff_operation( iff_bulk_list_t* bulk_list, char* iff, int line_no)
{
	iff_bulk_t* bulk1 = NULL;
	iff_bulk_t* bulk2 = NULL;
	char *pointer;

	if( iff_is_parent( iff )){
		bulk1 = iff_bulk_create();
		iff_bulk_regist( bulk1, iff, line_no );
		iff_bulk_list_regist( bulk_list, bulk1);
	} else if( iff_is_continuous_line( iff )){
		bulk1 = iff_bulk_list_search_pointing_bulk( bulk_list, iff_get_continuous_line_pointer(iff));
		if( bulk1 ) {
			iff_bulk_regist( bulk1, iff, line_no );
		} else {
			bulk1 = iff_bulk_create();
			iff_bulk_regist( bulk1, iff, line_no );
			iff_bulk_list_regist( bulk_list, bulk1);
		}
	} else {
		HECMW_assert(0);
	}


	if( iff_is_pointing_line( iff )){
		pointer = iff_get_pointing_line_pointer(iff);
		bulk2 = iff_bulk_list_search_continuous_bulk( bulk_list, pointer);
		if( bulk2 ){
			iff_bulk_move_append( bulk1, bulk2 );
			iff_bulk_list_remove( bulk_list, bulk2 );
			iff_bulk_free( bulk2 );
		}
	}

	if( iff_bulk_is_completed( bulk1 )){
		int rcode = iff_bulk_parse( bulk1 );
		iff_bulk_list_remove( bulk_list, bulk1 );
		iff_bulk_free( bulk1 );
		return rcode;
	}

	return 0;
}

/* f_open_close.c */ 

/* JP-444 */


/* JP-445 */

/* JP-446 */

static void file_stack_init()
{
	file_stack_pos = 0;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-447 */

static int file_stack_push(char* fname, FILE* fp, int lineno)
{
	if( file_stack_pos >= FILE_STACK_SIZE){
		/* JP-448 */
		return -1;
	}

	file_stack[ file_stack_pos ].fp = fp;
	strcpy( file_stack[ file_stack_pos ].filename, fname );
	file_stack[ file_stack_pos ].lineno = lineno;

	file_stack_pos++;
	return 0;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-449 */

static int file_stack_pop(char* fname, FILE** fp, int* lineno)
{
	if( file_stack_pos <= 0){
		/* JP-450 */
		return -1;
	}

	file_stack_pos--;

	*fp = file_stack[ file_stack_pos ].fp;
	strcpy( fname, file_stack[ file_stack_pos ].filename );
	*lineno = file_stack[ file_stack_pos ].lineno;

	return 0;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-451 */

static int file_stack_clear()
{
	int i;
	int fg;

	for(i = 0; i<file_stack_pos; i++){
		fclose( file_stack[i].fp );
	}

	fg =  (file_stack_pos != 0);
	file_stack_pos = 0;

	return fg;
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-452 */

static int file_stack_check( const char* fname )
{
	int i;
	for(i = 0; i<file_stack_pos; i++){
		if(strcmp( fname, file_stack[i].filename ) == 0)
			return -1;
	}
	return 0;
}



/* JP-453 */
/* JP-454 */





/* ------------------------------------------------------------------------------------------------ */
/* JP-455 */

static void nastran_file_init()
{
	c_fp = 0;
	c_filename[0] = 0;
	c_lineno = 0;

	file_stack_init();
}


/* ------------------------------------------------------------------------------------------------ */
/* JP-456 */

static int nastran_file_open( const char* filename )
{
	FILE* fp;

	if(file_stack_check( filename )){
		/* JP-457 */
		printf("エラー：多重オープンをしようとした\n");
		return -1;
	}

	fp = fopen( filename, "r");

	if(!fp){
		/* JP-458 */
		/* HECMW_set_error(HECMW_IO_NASTRAN_XXXX, "File: %s, %s", filename, strerror(errno)); */
		return -1;
	}

	if( c_fp ){
		if( file_stack_push( c_filename, c_fp, c_lineno )) { /* JP-459 */
			/* JP-460 */
			/* HECMW_set_error(HECMW_IO_NASTRAN_XXXX, "File: %s, %s", filename, strerror(errno)); */
			return -1;
		}
	}

	c_fp = fp;
	strcpy( c_filename, filename );
	c_lineno = 0;

	HECMW_io_set_gridfile((char*)filename);

	return 0;
}



/* ------------------------------------------------------------------------------------------------ */
/* JP-461 */


static int nastran_file_close()
{
	if(!c_fp){
		/* JP-462 */
		return -1;
	}

	if(fclose(c_fp)) {
		/* JP-463 */
		/* HECMW_set_error(HECMW_IO_HEC_E0002, "File: %s, %s", filename, strerror(errno)); */
		return -1;
	}

	if(file_stack_pop( c_filename, &c_fp, &c_lineno )) {
		HECMW_io_set_gridfile( "Unknown" );
		c_fp = NULL;
		strcpy( c_filename, "Unknown");
		return 1;
	} else {
		return 0;
	}
}







/* read_nastran.c */ 

/* JP-464 */
/* JP-465 */

/* ------------------------------------------------------------------------------------------------ */
/* JP-466 */

/* JP-467 */

static
int read_nastran( const char* filename )
{
	/* JP-468 */ 

	int end_fg;
	int read_fg;
	int next_read_fg;

	char iff_buff1[IFF_SIZE];
	char iff_buff2[IFF_SIZE];

	char* iff1;
	char* iff2;
	char* temp;

	int iff1_lineno;

	iff_bulk_list_t bulk_list;
	int bulk_number;

	int auto_pointer_counter = 1;

	iff_bulk_list_init( &bulk_list );

	if( nastran_file_open( filename ))
		return -1;

	skip_to_begin_bulk(c_fp, &c_lineno);		/* JP-469 */

	iff1 = iff_buff1;
	iff2 = iff_buff2;

	read_fg = read_iff( c_fp, iff1, &c_lineno);

	iff1_lineno = c_lineno;

	end_fg = 1;

	while(end_fg){
		next_read_fg = read_iff( c_fp, iff2, &c_lineno);
		switch( next_read_fg ) {
		case read_iff_success:
			if( iff_add_auto_pointer( iff1, iff2, &auto_pointer_counter ) )
					return -1;
			break;
		case read_iff_comment:
			continue;
		}
		switch( read_fg  ){
		case read_iff_success:
			iff_operation( &bulk_list, iff1, iff1_lineno ); /* JP-470 */
			break;
		case read_iff_error:
			return -1;
			break;
		case read_iff_comment:
			break;
		case read_iff_eof:
			if( nastran_file_close() != 0)
				end_fg = 0;
		}

		temp = iff1;
		iff1 = iff2;
		iff2 = temp;
		iff1_lineno = c_lineno;
		read_fg = next_read_fg;
	}

	iff_bulk_list_after_operation( &bulk_list );

	iff_bulk_list_clear( &bulk_list );

	return 0;
}
/* read_file.c */ 

/* JP-471 */





int HECMW_read_nastran_mesh(const char *filename)
{
	/* -----------------------------------

	HECMW_log(HECMW_LOG_DEBUG, "Start to read NASTRAN mesh");

	if(filename == NULL) { 
		HECMW_set_error(HECMW_IO_E0001, "Not specified filename for NASTRAN mesh input routine");
		return -1;
	}
	HECMW_log(HECMW_LOG_DEBUG, "NASTRAN mesh file is '%s'", filename);

	if(strlen(filename) > HECMW_FILENAME_LEN) {
		HECMW_set_error(HECMW_IO_E0002, "");
		return -1;
	}

	strcpy(grid_filename, filename);
	HECMW_io_set_gridfile(grid_filename);

	---------------------------------------- */
	
	nastran_file_init();
	section_list_init();

	if( read_nastran( filename )) return -1;

	if( surface_elem_type_decide())	return -1;
	if( section_list_finalize()) return -1;

	HECMW_io_set_gridfile((char*)filename);
	
	return 0;
}



/* ------------------------------------------------------------------------------------------------ */

struct hecmwST_local_mesh *
HECMW_get_nastran_mesh(const char *filename)
{
	struct hecmwST_local_mesh *local_mesh;

	if(HECMW_io_init()) return NULL;
	if(HECMW_io_pre_process()) return NULL;
	if(HECMW_read_nastran_mesh(filename)) return NULL;
	if(HECMW_io_post_process()) return NULL;
	local_mesh = HECMW_io_make_local_mesh();
	if(local_mesh == NULL) return NULL;
	if(HECMW_io_finalize()) return NULL;

	strcpy(grid_filename, "Unknown");

	return local_mesh;
}



