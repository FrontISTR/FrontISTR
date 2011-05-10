%{
#include <string.h>
#include <errno.h>
#include "hecmw_ablex.h"
#include "hecmw_util.h"

typedef union {
	double dval;
	char *str;
} YYSTYPE;

static YYSTYPE yylval;
static int lineno;
static int lineno_inc;
static int flag_including = 0;
static FILE *incfp;
static char include_filename[HECMW_FILENAME_LEN+1];
static YY_BUFFER_STATE prev_state;
static int flag_header;

static void set_flag_header(int flag);
%}

%option nounput 

comment		^(\*\*).*\r?\n
ws			[ \t]+
nl			\r?\n

int 		[+-]?[0-9]+
double		[+-]?[0-9]*\.([0-9]*([eE][+-]?[0-9]+)?)?
name		[a-zA-Z][_a-zA-Z0-9-]*
file		([a-zA-Z]:)*[\\/._a-zA-Z0-9-]+

%x HEADER
%s INCLUDE

%%

<HEADER>{comment}	{ lineno++; set_flag_header(0); }
<HEADER>{nl}		{ lineno++; set_flag_header(0); }
<HEADER>"*HEADING".*	{ return flag_header ? HECMW_ABLEX_H_HEADING : ' '; }
<HEADER>^"*".*		{
						HECMW_ablex_unput_token();
						set_flag_header(1);
						BEGIN INITIAL;
					}
<HEADER>^.*			{ set_flag_header(0); return HECMW_ABLEX_HEADER; }

<INCLUDE>{nl}		{ lineno_inc++; set_flag_header(0); return HECMW_ABLEX_NL; }
<INCLUDE><<EOF>>	{
						yy_delete_buffer(YY_CURRENT_BUFFER);
						yy_switch_to_buffer(prev_state);
						if(fclose(incfp)) {
							HECMW_set_error(HECMW_IO_E0006,
								"File: %s, %s", include_filename, strerror(errno));
						}
						flag_including = 0;
						BEGIN INITIAL;
					}

{comment}	{ lineno++; set_flag_header(0); }
{ws}		{ set_flag_header(0); }
{nl}		{ lineno++; set_flag_header(0); return HECMW_ABLEX_NL; }

^"*"AMPLITUDE	{ set_flag_header(1); return HECMW_ABLEX_H_AMPLITUDE; }
"*"AMPLITUDE	{ return flag_header ? HECMW_ABLEX_H_AMPLITUDE : ' '; }
^"*"CONDUCTIVITY	{ set_flag_header(1); return HECMW_ABLEX_H_CONDUCTIVITY; }
"*"CONDUCTIVITY		{ return flag_header ? HECMW_ABLEX_H_CONDUCTIVITY : ' '; }
^"*"DENSITY	{ set_flag_header(1); return HECMW_ABLEX_H_DENSITY; }
"*"DENSITY	{ return flag_header ? HECMW_ABLEX_H_DENSITY : ' '; }
^"*"ELASTIC	{ set_flag_header(1); return HECMW_ABLEX_H_ELASTIC; }
"*"ELASTIC	{ return flag_header ? HECMW_ABLEX_H_ELASTIC : ' '; }
^"*"ELEMENT	{ set_flag_header(1); return HECMW_ABLEX_H_ELEMENT; }
"*"ELEMENT	{ return flag_header ? HECMW_ABLEX_H_ELEMENT : ' '; }
^"*"ELSET	{ set_flag_header(1); return HECMW_ABLEX_H_ELSET; }
"*"ELSET	{ return flag_header ? HECMW_ABLEX_H_ELSET : ' '; }
^"*"EQUATION	{ set_flag_header(1); return HECMW_ABLEX_H_EQUATION; }
"*"EQUATION	{ return flag_header ? HECMW_ABLEX_H_EQUATION : ' '; }
^"*"HEADING	{
				BEGIN HEADER; 
				set_flag_header(1);
				return HECMW_ABLEX_H_HEADING;
			}
^"*"INCLUDE	{ set_flag_header(1); return HECMW_ABLEX_H_INCLUDE; }
"*"INCLUDE	{ return flag_header ? HECMW_ABLEX_H_INCLUDE : ' '; }
^"*"INITIAL" "?CONDITIONS	{ set_flag_header(1); return HECMW_ABLEX_H_INITIAL; }
"*"INITIAL" "?CONDITIONS	{ return flag_header ? HECMW_ABLEX_H_INITIAL : ' '; }
^"*"MATERIAL	{ set_flag_header(1); return HECMW_ABLEX_H_MATERIAL; }
"*"MATERIAL		{ return flag_header ? HECMW_ABLEX_H_MATERIAL : ' '; }
^"*"NODE	{ set_flag_header(1); return HECMW_ABLEX_H_NODE; }
"*"NODE		{ return flag_header ? HECMW_ABLEX_H_NODE : ' '; }
^"*"NSET	{ set_flag_header(1); return HECMW_ABLEX_H_NSET; }
"*"NSET		{ return flag_header ? HECMW_ABLEX_H_NSET : ' '; }
^"*"SHELL" "?SECTION	{ set_flag_header(1); return HECMW_ABLEX_H_SHELL_SECTION; }
"*"SHELL" "?SECTION		{ return flag_header ? HECMW_ABLEX_H_SHELL_SECTION : ' '; }
^"*"SOLID" "?SECTION	{ set_flag_header(1); return HECMW_ABLEX_H_SOLID_SECTION; }
"*"SOLID" "?SECTION		{ return flag_header ? HECMW_ABLEX_H_SOLID_SECTION : ' '; }
^"*"SPECIFIC" "?HEAT	{ set_flag_header(1); return HECMW_ABLEX_H_SPECIFIC_HEAT; }
"*"SPECIFIC" "?HEAT		{ return flag_header ? HECMW_ABLEX_H_SPECIFIC_HEAT : ' '; }
^"*"SYSTEM	{ set_flag_header(1); return HECMW_ABLEX_H_SYSTEM; }
"*"SYSTEM	{ return flag_header ? HECMW_ABLEX_H_SYSTEM : ' '; }

ABSOLUTE	{ set_flag_header(0); return HECMW_ABLEX_K_ABSOLUTE; }
ANISOTROPIC	{ set_flag_header(0); return HECMW_ABLEX_K_ANISOTROPIC; }
DEFINITION	{ set_flag_header(0); return HECMW_ABLEX_K_DEFINITION; }
ELSET		{ set_flag_header(0); return HECMW_ABLEX_K_ELSET; }
ENGINEERING" "?CONSTANTS	{ set_flag_header(0); return HECMW_ABLEX_K_ENGINEERING_CONSTANTS; }
GENERATE	{ set_flag_header(0); return HECMW_ABLEX_K_GENERATE; }
INPUT		{ set_flag_header(0); return HECMW_ABLEX_K_INPUT; }
ISOTROPIC	{ set_flag_header(0); return HECMW_ABLEX_K_ISOTROPIC; }
LAMINA		{ set_flag_header(0); return HECMW_ABLEX_K_LAMINA; }
MATERIAL	{ set_flag_header(0); return HECMW_ABLEX_K_MATERIAL; }
NAME		{ set_flag_header(0); return HECMW_ABLEX_K_NAME; }
NSET		{ set_flag_header(0); return HECMW_ABLEX_K_NSET; }
ORTHOTROPIC	{ set_flag_header(0); return HECMW_ABLEX_K_ORTHOTROPIC; }
RELATIVE	{ set_flag_header(0); return HECMW_ABLEX_K_RELATIVE; }
STEP" "?TIME	{ set_flag_header(0); return HECMW_ABLEX_K_STEP_TIME; }
SYSTEM		{ set_flag_header(0); return HECMW_ABLEX_K_SYSTEM; }
TABULAR		{ set_flag_header(0); return HECMW_ABLEX_K_TABULAR; }
TEMPERATURE	{ set_flag_header(0); return HECMW_ABLEX_K_TEMPERATURE; }
TIME		{ set_flag_header(0); return HECMW_ABLEX_K_TIME; }
TYPE		{ set_flag_header(0); return HECMW_ABLEX_K_TYPE; }
UNSORTED	{ set_flag_header(0); return HECMW_ABLEX_K_UNSORTED; }
VALUE		{ set_flag_header(0); return HECMW_ABLEX_K_VALUE; }

B31		{ set_flag_header(0); return HECMW_ABLEX_E_B31; }
B32		{ set_flag_header(0); return HECMW_ABLEX_E_B32; }
C3D4	{ set_flag_header(0); return HECMW_ABLEX_E_C3D4; }
C3D6	{ set_flag_header(0); return HECMW_ABLEX_E_C3D6; }
C3D8	{ set_flag_header(0); return HECMW_ABLEX_E_C3D8; }
C3D8I	{ set_flag_header(0); return HECMW_ABLEX_E_C3D8I; }
C3D10	{ set_flag_header(0); return HECMW_ABLEX_E_C3D10; }
C3D15	{ set_flag_header(0); return HECMW_ABLEX_E_C3D15; }
C3D20	{ set_flag_header(0); return HECMW_ABLEX_E_C3D20; }
CAX3	{ set_flag_header(0); return HECMW_ABLEX_E_CAX3; }
CAX4	{ set_flag_header(0); return HECMW_ABLEX_E_CAX4; }
CAX4I	{ set_flag_header(0); return HECMW_ABLEX_E_CAX4I; }
CAX4R	{ set_flag_header(0); return HECMW_ABLEX_E_CAX4R; }
CAX6	{ set_flag_header(0); return HECMW_ABLEX_E_CAX6; }
CAX8	{ set_flag_header(0); return HECMW_ABLEX_E_CAX8; }
CAX8R	{ set_flag_header(0); return HECMW_ABLEX_E_CAX8R; }
CPE3	{ set_flag_header(0); return HECMW_ABLEX_E_CPE3; }
CPE4	{ set_flag_header(0); return HECMW_ABLEX_E_CPE4; }
CPE4I	{ set_flag_header(0); return HECMW_ABLEX_E_CPE4I; }
CPE4R	{ set_flag_header(0); return HECMW_ABLEX_E_CPE4R; }
CPE6	{ set_flag_header(0); return HECMW_ABLEX_E_CPE6; }
CPE8	{ set_flag_header(0); return HECMW_ABLEX_E_CPE8; }
CPE8R	{ set_flag_header(0); return HECMW_ABLEX_E_CPE8R; }
CPS3	{ set_flag_header(0); return HECMW_ABLEX_E_CPS3; }
CPS4	{ set_flag_header(0); return HECMW_ABLEX_E_CPS4; }
CPS4I	{ set_flag_header(0); return HECMW_ABLEX_E_CPS4I; }
CPS4R	{ set_flag_header(0); return HECMW_ABLEX_E_CPS4R; }
CPS6	{ set_flag_header(0); return HECMW_ABLEX_E_CPS6; }
CPS8	{ set_flag_header(0); return HECMW_ABLEX_E_CPS8; }
CPS8R	{ set_flag_header(0); return HECMW_ABLEX_E_CPS8R; }
DC1D2	{ set_flag_header(0); return HECMW_ABLEX_E_DC1D2; }
DC1D3	{ set_flag_header(0); return HECMW_ABLEX_E_DC1D3; }
DC2D3	{ set_flag_header(0); return HECMW_ABLEX_E_DC2D3; }
DC2D4	{ set_flag_header(0); return HECMW_ABLEX_E_DC2D4; }
DC2D6	{ set_flag_header(0); return HECMW_ABLEX_E_DC2D6; }
DC2D8	{ set_flag_header(0); return HECMW_ABLEX_E_DC2D8; }
DC3D4	{ set_flag_header(0); return HECMW_ABLEX_E_DC3D4; }
DC3D6	{ set_flag_header(0); return HECMW_ABLEX_E_DC3D6; }
DC3D8	{ set_flag_header(0); return HECMW_ABLEX_E_DC3D8; }
DC3D10	{ set_flag_header(0); return HECMW_ABLEX_E_DC3D10; }
DC3D15	{ set_flag_header(0); return HECMW_ABLEX_E_DC3D15; }
DC3D20	{ set_flag_header(0); return HECMW_ABLEX_E_DC3D20; }
DCAX3	{ set_flag_header(0); return HECMW_ABLEX_E_DCAX3; }
DCAX4	{ set_flag_header(0); return HECMW_ABLEX_E_DCAX4; }
DCAX6	{ set_flag_header(0); return HECMW_ABLEX_E_DCAX6; }
DCAX8	{ set_flag_header(0); return HECMW_ABLEX_E_DCAX8; }
DINTER4	{ set_flag_header(0); return HECMW_ABLEX_E_DINTER4; }
DINTER8	{ set_flag_header(0); return HECMW_ABLEX_E_DINTER8; }
DS4		{ set_flag_header(0); return HECMW_ABLEX_E_DS4; }
DS8		{ set_flag_header(0); return HECMW_ABLEX_E_DS8; }
INTER4	{ set_flag_header(0); return HECMW_ABLEX_E_INTER4; }
INTER8	{ set_flag_header(0); return HECMW_ABLEX_E_INTER8; }
S3R		{ set_flag_header(0); return HECMW_ABLEX_E_S3R; }
S4R		{ set_flag_header(0); return HECMW_ABLEX_E_S4R; }
S8R		{ set_flag_header(0); return HECMW_ABLEX_E_S8R; }
T3D2	{ set_flag_header(0); return HECMW_ABLEX_E_T3D2; }
T3D3	{ set_flag_header(0); return HECMW_ABLEX_E_T3D3; }

=			{ set_flag_header(0); return '='; }
,			{ set_flag_header(0); return ','; }
C			{ set_flag_header(0); return 'C'; }
R			{ set_flag_header(0); return 'R'; }

{int}		{
				yylval.dval = atof(yytext);
				set_flag_header(0);
				return HECMW_ABLEX_INT;
			}
{double}	{
				yylval.dval = atof(yytext);
				set_flag_header(0);
				return HECMW_ABLEX_DOUBLE;
			}
{name}		{
				yylval.str = yytext;
				set_flag_header(0);
				return HECMW_ABLEX_NAME;
			}
{file}		{
				yylval.str = yytext;
				set_flag_header(0);
				return HECMW_ABLEX_FILENAME;
			}
.			{ set_flag_header(0); return yytext[0]; }
<<EOF>>		{ return 0; }

%%

static void
set_flag_header(int flag)
{
	flag_header = flag ? 1 : 0;
}


int
HECMW_ablex_get_lineno(void)
{
	return flag_including ? lineno_inc : lineno;
}


double
HECMW_ablex_get_number(void)
{
	return yylval.dval;
}


char *
HECMW_ablex_get_text(void)
{
	return yytext;
}


int
HECMW_ablex_is_including(void)
{
	return flag_including ? 1 : 0;
}


int
HECMW_ablex_next_token(void)
{
	return yylex();
}


int
HECMW_ablex_next_token_skip(int skip_token)
{
	int token;
	while((token = yylex())) {
		if(token != skip_token) break;
	}
	return token;
}


int
HECMW_ablex_set_input(FILE *fp)
{
	static int first = 1;	
	if(fp == NULL) return -1;
	if(first) {
		yyin = fp;
		first = 0;
	} else {
		yyrestart(fp);
	}
	lineno = 1;
	return 0;
}


int
HECMW_ablex_skip_line(void)
{
	int token;
	while((token = HECMW_ablex_next_token())) {
		if(token == HECMW_ABLEX_NL) break;
	}
	return token;
}


int
HECMW_ablex_switch_to_include(const char *filename)
{
	if(filename == NULL) return -1;
	if(flag_including) {
		HECMW_set_error(HECMW_IO_E0004, "File: %s", filename);
		return -1;
	}
	if((incfp = fopen(filename, "r")) == NULL) {
		HECMW_set_error(HECMW_IO_E0005, "File: %s, %s", filename, strerror(errno));
		return -1;
	}
	strcpy(include_filename, filename);
	prev_state = YY_CURRENT_BUFFER;
	yy_switch_to_buffer(yy_create_buffer(incfp, YY_BUF_SIZE));
	lineno_inc = 1;
	flag_including = 1;
	BEGIN INCLUDE;
	return 0;
}


int
HECMW_ablex_unput_token(void)
{
	char *p = yytext;
	while(*p) {
		if(*p == '\n') {
			if(YYSTATE == INCLUDE) {
				lineno_inc--;
			} else {
				lineno--;
			}
		}
		p++;
	}
	yyless(0);
	return 0;
}


int
HECMW_abwrap(void)
{
	return 1;
}


