/*=====================================================================*
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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

/**
 * @file fstr_ctrl_util.h
 * Utility for open/close and parse of fstr control file
 * caution: define parameter 'ctrl' as integer for call from fortran
 * @date 2004/08/11
 * @author Noboru Imai
*/


#ifndef fstr_ctrl_utilH
#define fstr_ctrl_utilH

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <hecmw_config.h>


#define FSTR_CTRL_RCODE_PARAM_SUCCESS       0
#define FSTR_CTRL_RCODE_PARAM_ERROR        -1
#define FSTR_CTRL_RCODE_PARAM_TYPE_ERROR   -2
#define FSTR_CTRL_RCODE_PARAM_RANGE_ERROR  -3
#define FSTR_CTRL_RCODE_PARAM_NOTHING       1
#define FSTR_CTRL_RCODE_PARAM_VALUE_NOTHING 2

#define FSTR_CTRL_RCODE_DATA_SUCCESS        0
#define FSTR_CTRL_RCODE_DATA_ERROR         -1
#define FSTR_CTRL_RCODE_DATA_TYPE_ERROR    -2
#define FSTR_CTRL_RCODE_DATA_RANGE_ERROR   -3
#define FSTR_CTRL_RCODE_DATA_NOTHING        1
#define FSTR_CTRL_RCODE_DATA_LINE_NOTHING   2


typedef struct st_ctrl_rec {
	int line_no;
	char* line;
} ctrl_rec;


typedef struct st_fstr_ctrl_data {
	ctrl_rec* rec;
	int rec_n;
	int* header_pos;
	int header_n;
	int* data_line_n;
	int current_header_index;
} fstr_ctrl_data;


/* ctrl list for fortran interface  */
#define ctrl_list_size 20
#ifdef fstr_ctrl_util_MAIN
fstr_ctrl_data* ctrl_list[ ctrl_list_size ];
#else
extern fstr_ctrl_data* ctrl_list[ ctrl_list_size ];
#endif


/* ================================================================================= */
/**
 * Get error message (for _ex function )
 * @param buff buffer for copied message
 */

void fstr_ctrl_get_err_msg(char* buff);


/* ================================================================================= */
/**
 * Open FSTR control file
 * @param filename Specify FSTR control file name (included path)
 * @return fstr_ctrl_data* pointer or NULL
 */

fstr_ctrl_data* fstr_ctrl_open( const char* filename );


/* ================================================================================= */
/**
 * Obtaining record number
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @return record number without comment and blank lines
 */

int fstr_ctrl_get_rec_number( fstr_ctrl_data* ctrl );

/* ================================================================================= */
/**
 * Obtaining line
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @param rec_no Specify record no
 * @param buff Pointer to copy the record string line
 * @return record number without comment and blank lines
 */

int fstr_ctrl_get_line( fstr_ctrl_data* ctrl, int rec_no, char* buff );

/* ================================================================================= */
/**
 * Seeking to specified header (current header)
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @param header_name Specify header name or NULL pointer to seek next header
 *        When this function is called from fortran, use blank string instead of NULL
 * @return 1:success, 0:fail ( REMARK! )
 */

int fstr_ctrl_seek_header( fstr_ctrl_data* ctrl, const char* header_name );

/* ================================================================================= */
/**
 * Seeking next header (current header)
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @return 1:success, 0:fail ( REMARK! )
 */

int fstr_ctrl_seek_next_header( fstr_ctrl_data* ctrl );

/* ================================================================================= */
/**
 * Obtaining current header name
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @param header_name buffer for header name
 * @return 0:success, -1:fail
 */

int fstr_ctrl_get_current_header_name( fstr_ctrl_data* ctrl, char* header_name );


/* ================================================================================= */
/**
 * Obtaining line number of current header in fstr control file
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @return line number or 0:not header seeking, -1:error
*/

int fstr_ctrl_get_current_header_line_no( fstr_ctrl_data* ctrl );

/* ================================================================================= */
/**
 * Obtaining record line number of current header in fstr control file
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @return line number or 0:not header seeking, -1:error
*/

int fstr_ctrl_get_current_header_pos( fstr_ctrl_data* ctrl );

/* ================================================================================= */
/**
 * Obtaining value of parameter in current header line
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @param paran_name Specify parameter name
 * @param value_list Specify value list by csv string
 * @param type Specify type of parameter's value as charactor
 *        'I':integer(int), 'C' or 'S':string(char*), 'R':real(double),
 *        'P':Pattern Input(first string in value_list is 0) (int)
 *        'E':exist or not(int, 1 or 0)
 * @param val Put pointer to store the parameter's value
 * @return 0:success, 1:nothing parameter, 2:nothing value of param., 3:type change error
 *         4:value range error
 */


int fstr_ctrl_get_param( fstr_ctrl_data* ctrl, const char* param_name, const char* value_list, char type, void* val );


/* ================================================================================= */
/**
 * Obtaining value of parameter in current header line
 * This function uses "fstr_ctrl_get_param" and prints message in error occurance.
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @param paran_name Specify parameter name
 * @param necessity Necessity of parameter (1 or 0)
 * @param type Specify type of parameter's value as charactor
 *        'I':integer(int), 'C' or 'S':string(char*), 'R':real(double),
 *        'P':Pattern Input(first string in value_list is 1) (int)
 *        'E':exist or not(int, 1 or 0)
 * @param val Put pointer to store the parameter's value
 * @return 0:success, 1:nothing parameter, 2:nothing value of param., 3:type change error
 *         4:value range error
 */

int fstr_ctrl_get_param_ex( fstr_ctrl_data* ctrl,
	const char* param_name, const char* value_list, int necessity, char type, void* val );

/* ================================================================================= */
/**
 * Obtaining data line number of current header
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @return Data line number or -1 (Error)
 */

int fstr_ctrl_get_data_line_n( fstr_ctrl_data* ctrl );

/* ================================================================================= */
/* JP-0 */

int fstr_ctrl_copy_data_line( fstr_ctrl_data* ctrl, int line_no, char* data_line );

/* ================================================================================= */
/* JP-1 */

int fstr_ctrl_get_data_n_in_line( fstr_ctrl_data* ctrl, int line_no, const char* delim );


/* ================================================================================= */
/**
 * Obtainig converting error position in previous fstr_ctrl_get_data executation
 * @return Position of error, 0: no error, or, -1 : non converting error
 */

int fstr_ctrl_get_data_error_pos(void);

int fstr_ctrl_get_data_error_line(void);

/* ================================================================================= */
/* JP-2 */

int fstr_ctrl_get_data( fstr_ctrl_data* ctrl, int line_no, const char* format, ... );
int fstr_ctrl_get_data_v( fstr_ctrl_data* ctrl, int line_no, const char* format, va_list va );
int fstr_ctrl_get_data_ex( fstr_ctrl_data* ctrl, int line_no, const char* format, ... );

int fstr_ctrl_get_data_array( fstr_ctrl_data* ctrl, const char* format, ... );
int fstr_ctrl_get_data_array_v( fstr_ctrl_data* ctrl, const char* format, va_list va );
int fstr_ctrl_get_data_array_ex( fstr_ctrl_data* ctrl, const char* format, ... );


/* ================================================================================= */
/**
 * Closing fstr control file
 * @param ctrl Specify fstr_ctrl_data* pointer obtained from fstr_ctrl_open
 * @return 0:success, -1:fail
 */

int fstr_ctrl_close( fstr_ctrl_data* ctrl );


#endif
