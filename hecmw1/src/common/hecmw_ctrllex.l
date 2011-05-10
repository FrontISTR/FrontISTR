%{
#include <string.h>
#include <errno.h>
#include "hecmw_ctrllex.h"
#include "hecmw_util.h"

typedef union {
	double dval;
	char *str;
} YYSTYPE;

static YYSTYPE yylval;
static int lineno;
static int flag_header;

static void set_flag_header(int flag);
%}

%option nounput

comment		^(!!|#).*\r?\n
ws			[ \t]+
nl			\r?\n

int 		[+-]?[0-9]+
double		[+-]?[0-9]*\.([0-9]*([eE][+-]?[0-9]+)?)?
name		[_a-zA-Z][_a-zA-Z0-9-]*
file		([a-zA-Z]:)*[\\/._a-zA-Z0-9-]+

%%

{comment}	{ lineno++; set_flag_header(0); }
{ws}		{ set_flag_header(0); }
{nl}		{ lineno++; set_flag_header(0); return HECMW_CTRLLEX_NL; }

^!CONTROL	{ set_flag_header(1); return HECMW_CTRLLEX_H_CONTROL; }
!CONTROL	{ return flag_header ? HECMW_CTRLLEX_H_CONTROL : ' '; }
^!MESH		{ set_flag_header(1); return HECMW_CTRLLEX_H_MESH; }
!MESH		{ return flag_header ? HECMW_CTRLLEX_H_MESH : ' '; }
^!MESH" "?GROUP	{ set_flag_header(1); return HECMW_CTRLLEX_H_MESH_GROUP; }
!MESH" "?GROUP	{ return flag_header ? HECMW_CTRLLEX_H_MESH_GROUP : ' '; }
^!RESULT	{ set_flag_header(1); return HECMW_CTRLLEX_H_RESULT; }
!RESULT		{ return flag_header ? HECMW_CTRLLEX_H_RESULT : ' '; }
^!RESTART	{ set_flag_header(1); return HECMW_CTRLLEX_H_RESTART; }
!RESTART	{ return flag_header ? HECMW_CTRLLEX_H_RESTART : ' '; }
^!COUPLE" "?UNIT	{ set_flag_header(1); return HECMW_CTRLLEX_H_COUPLE_UNIT; }
!COUPLE" "?UNIT		{ return flag_header ? HECMW_CTRLLEX_H_COUPLE_UNIT : ' '; }
^!COUPLE	{ set_flag_header(1); return HECMW_CTRLLEX_H_COUPLE; }
!COUPLE	{ return flag_header ? HECMW_CTRLLEX_H_COUPLE : ' '; }
^!COUPLE" "?BOUNDARY	{ set_flag_header(1); return HECMW_CTRLLEX_H_COUPLE_BOUNDARY; }
!COUPLE" "?BOUNDARY		{ return flag_header ? HECMW_CTRLLEX_H_COUPLE_BOUNDARY : ' '; }
^!UNIT1	{ set_flag_header(1); return HECMW_CTRLLEX_H_UNIT1; }
!UNIT1	{ return flag_header ? HECMW_CTRLLEX_H_UNIT1 : ' '; }
^!UNIT2	{ set_flag_header(1); return HECMW_CTRLLEX_H_UNIT2; }
!UNIT2	{ return flag_header ? HECMW_CTRLLEX_H_UNIT2 : ' '; }

ABAQUS			{ set_flag_header(0); return HECMW_CTRLLEX_K_ABAQUS; }
DIR				{ set_flag_header(0); return HECMW_CTRLLEX_K_DIR; }
FEMAP			{ set_flag_header(0); return HECMW_CTRLLEX_K_FEMAP; }
GeoFEM			{ set_flag_header(0); return HECMW_CTRLLEX_K_GEOFEM; }
HECMW-DIST		{ set_flag_header(0); return HECMW_CTRLLEX_K_HECMW_DIST; }
HECMW-ENTIRE	{ set_flag_header(0); return HECMW_CTRLLEX_K_HECMW_ENTIRE; }
IN				{ set_flag_header(0); return HECMW_CTRLLEX_K_IN; }
INOUT			{ set_flag_header(0); return HECMW_CTRLLEX_K_INOUT; }
IO				{ set_flag_header(0); return HECMW_CTRLLEX_K_IO; }
REFINE			{ set_flag_header(0); return HECMW_CTRLLEX_K_REFINE; }
NAME			{ set_flag_header(0); return HECMW_CTRLLEX_K_NAME; }
NASTRAN			{ set_flag_header(0); return HECMW_CTRLLEX_K_NASTRAN; }
OUT				{ set_flag_header(0); return HECMW_CTRLLEX_K_OUT; }
TYPE			{ set_flag_header(0); return HECMW_CTRLLEX_K_TYPE; }
NPROC			{ set_flag_header(0); return HECMW_CTRLLEX_K_NPROC; }
MXN				{ set_flag_header(0); return HECMW_CTRLLEX_K_MXN; }
MAXMN			{ set_flag_header(0); return HECMW_CTRLLEX_K_MAXMN; }
MANUAL			{ set_flag_header(0); return HECMW_CTRLLEX_K_MANUAL; }
UNIT1			{ set_flag_header(0); return HECMW_CTRLLEX_K_UNIT1; }
UNIT2			{ set_flag_header(0); return HECMW_CTRLLEX_K_UNIT2; }
COUPLE			{ set_flag_header(0); return HECMW_CTRLLEX_K_COUPLE; }
DIRECTION		{ set_flag_header(0); return HECMW_CTRLLEX_K_DIRECTION; }
UNIT1_TO_UNIT2	{ set_flag_header(0); return HECMW_CTRLLEX_K_UNIT1_TO_UNIT2; }
UNIT2_TO_UNIT1	{ set_flag_header(0); return HECMW_CTRLLEX_K_UNIT2_TO_UNIT1; }
INTERPOLATION	{ set_flag_header(0); return HECMW_CTRLLEX_K_INTERPOLATION; }
TOLERANCE		{ set_flag_header(0); return HECMW_CTRLLEX_K_TOLERANCE; }
BBCOEF			{ set_flag_header(0); return HECMW_CTRLLEX_K_BBCOEF; }
BGCOEF			{ set_flag_header(0); return HECMW_CTRLLEX_K_BGCOEF; }
GEOM			{ set_flag_header(0); return HECMW_CTRLLEX_K_GEOM; }
DATA			{ set_flag_header(0); return HECMW_CTRLLEX_K_DATA; }
NODE			{ set_flag_header(0); return HECMW_CTRLLEX_K_NODE; }
SURFACE			{ set_flag_header(0); return HECMW_CTRLLEX_K_SURFACE; }
ELEMENT			{ set_flag_header(0); return HECMW_CTRLLEX_K_ELEMENT; }


=			{ set_flag_header(0); return '='; }
,			{ set_flag_header(0); return ','; }

{int}		{
				yylval.dval = atof(yytext);
				set_flag_header(0);
				return HECMW_CTRLLEX_INT;
			}
{double}	{
				yylval.dval = atof(yytext);
				set_flag_header(0);
				return HECMW_CTRLLEX_DOUBLE;
			}
{name}		{
				yylval.str = yytext;
				set_flag_header(0);
				return HECMW_CTRLLEX_NAME;
			}
{file}		{
				yylval.str = yytext;
				set_flag_header(0);
				return HECMW_CTRLLEX_FILENAME;
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
HECMW_ctrllex_get_lineno(void)
{
	return lineno;
}


double
HECMW_ctrllex_get_number(void)
{
	return yylval.dval;
}


char *
HECMW_ctrllex_get_text(void)
{
	return yytext;
}


int
HECMW_ctrllex_next_token(void)
{
	return yylex();
}


int
HECMW_ctrllex_next_token_skip(int skip_token)
{
	int token;
	while((token = yylex())) {
		if(token != skip_token) break;
	}
	return token;
}


int
HECMW_ctrllex_set_input(FILE *fp)
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
HECMW_ctrllex_skip_line(void)
{
	int token;
	while((token = HECMW_ctrllex_next_token())) {
		if(token == HECMW_CTRLLEX_NL) break;
	}
	return token;
}


int
HECMW_ctrllex_unput_token(void)
{
	char *p = yytext;
	while(*p) {
		if(*p == '\n') {
			lineno--;
		}
		p++;
	}
	yyless(0);
	return 0;
}


int
HECMW_ctrlwrap(void)
{
	return 1;
}


