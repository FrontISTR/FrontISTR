%{
#include <string.h>
#include <errno.h>
#include "hecmw_heclex.h"
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

comment		^(!!|#).*\r?\n
ws			[ \t]+
nl			\r?\n

int 		[+-]?[0-9]+
double		[+-]?([0-9]+\.?|\.?[0-9]+)([0-9]*([eE][+-]?[0-9]+)?)?
name		[_a-zA-Z][_a-zA-Z0-9-]*
file		([a-zA-Z]:)*[\\/._a-zA-Z0-9-]+

%x HEADER
%s INCLUDE

%%

<HEADER>{comment}	{ lineno++; set_flag_header(0); }
<HEADER>{nl}		{ lineno++; set_flag_header(0); }
<HEADER>!HEADER.*	{ return flag_header ? HECMW_HECLEX_H_HEADER : ' '; }
<HEADER>^!.*		{
						HECMW_heclex_unput_token();
						set_flag_header(1);
						BEGIN INITIAL;
					}
<HEADER>^.*			{ set_flag_header(0); return HECMW_HECLEX_HEADER; }

<INCLUDE>{nl}		{ lineno_inc++; set_flag_header(0); return HECMW_HECLEX_NL; }
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
{nl}		{ lineno++; set_flag_header(0); return HECMW_HECLEX_NL; }

^!AMPLITUDE	{ set_flag_header(1); return HECMW_HECLEX_H_AMPLITUDE; }
!AMPLITUDE	{ return flag_header ? HECMW_HECLEX_H_AMPLITUDE : ' '; }
^!CONNECTIVITY	{ set_flag_header(1); return HECMW_HECLEX_H_CONNECTIVITY; }
!CONNECTIVITY	{ return flag_header ? HECMW_HECLEX_H_CONNECTIVITY : ' '; }
^!CONTACT" "?PAIR	{ set_flag_header(1); return HECMW_HECLEX_H_CONTACT_PAIR; }
!CONTACT" "?PAIR	{ return flag_header ? HECMW_HECLEX_H_CONTACT_PAIR : ' '; }
^!EMBED" "?PAIR	{ set_flag_header(1); return HECMW_HECLEX_H_EMBED_PAIR; }
!EMBED" "?PAIR	{ return flag_header ? HECMW_HECLEX_H_EMBED_PAIR : ' '; }
^!ECOPY		{ set_flag_header(1); return HECMW_HECLEX_H_ECOPY; }
!ECOPY		{ return flag_header ? HECMW_HECLEX_H_ECOPY : ' '; }
^!EGEN		{ set_flag_header(1); return HECMW_HECLEX_H_EGEN; }
!EGEN		{ return flag_header ? HECMW_HECLEX_H_EGEN : ' '; }
^!EGROUP	{ set_flag_header(1); return HECMW_HECLEX_H_EGROUP; }
!EGROUP		{ return flag_header ? HECMW_HECLEX_H_EGROUP : ' '; }
^!ELEMENT	{ set_flag_header(1); return HECMW_HECLEX_H_ELEMENT; }
!ELEMENT	{ return flag_header ? HECMW_HECLEX_H_ELEMENT : ' '; }
^!END		{ set_flag_header(1); return HECMW_HECLEX_H_END; }
!END		{ return flag_header ? HECMW_HECLEX_H_END : ' '; }
^!EQUATION	{ set_flag_header(1); return HECMW_HECLEX_H_EQUATION; }
!EQUATION	{ return flag_header ? HECMW_HECLEX_H_EQUATION : ' '; }
^!HEADER	{
				BEGIN HEADER;
				set_flag_header(1);
				return HECMW_HECLEX_H_HEADER;
			}
^!INCLUDE	{ set_flag_header(1); return HECMW_HECLEX_H_INCLUDE; }
!INCLUDE	{ return flag_header ? HECMW_HECLEX_H_INCLUDE : ' '; }
^!INITIAL" "?CONDITION	{
							set_flag_header(1);
							return HECMW_HECLEX_H_INITIAL;
						}
!INITIAL" "?CONDITION	{
							return flag_header ? HECMW_HECLEX_H_INITIAL : ' ';
						}
^!ITEM		{ set_flag_header(1); return HECMW_HECLEX_H_ITEM; }
!ITEM		{ return flag_header ? HECMW_HECLEX_H_ITEM : ' '; }
^!MATERIAL	{ set_flag_header(1); return HECMW_HECLEX_H_MATERIAL; }
!MATERIAL	{ return flag_header ? HECMW_HECLEX_H_MATERIAL : ' '; }
^!NCOPY		{ set_flag_header(1); return HECMW_HECLEX_H_NCOPY; }
!NCOPY		{ return flag_header ? HECMW_HECLEX_H_NCOPY : ' '; }
^!NFILL		{ set_flag_header(1); return HECMW_HECLEX_H_NFILL; }
!NFILL		{ return flag_header ? HECMW_HECLEX_H_NFILL : ' '; }
^!NGEN		{ set_flag_header(1); return HECMW_HECLEX_H_NGEN; }
!NGEN		{ return flag_header ? HECMW_HECLEX_H_NGEN : ' '; }
^!NGROUP	{ set_flag_header(1); return HECMW_HECLEX_H_NGROUP; }
!NGROUP		{ return flag_header ? HECMW_HECLEX_H_NGROUP : ' '; }
^!NODE		{ set_flag_header(1); return HECMW_HECLEX_H_NODE; }
!NODE		{ return flag_header ? HECMW_HECLEX_H_NODE : ' '; }
^!SECTION	{ set_flag_header(1); return HECMW_HECLEX_H_SECTION; }
!SECTION	{ return flag_header ? HECMW_HECLEX_H_SECTION : ' '; }
^!SGROUP	{ set_flag_header(1); return HECMW_HECLEX_H_SGROUP; }
!SGROUP		{ return flag_header ? HECMW_HECLEX_H_SGROUP : ' '; }
^!SYSTEM	{ set_flag_header(1); return HECMW_HECLEX_H_SYSTEM; }
!SYSTEM		{ return flag_header ? HECMW_HECLEX_H_SYSTEM : ' '; }
^!ZERO		{ set_flag_header(1); return HECMW_HECLEX_H_ZERO; }
!ZERO		{ return flag_header ? HECMW_HECLEX_H_ZERO : ' '; }


ABAQUS		{ set_flag_header(0); return HECMW_HECLEX_K_ABAQUS; }
ABSOLUTE	{ set_flag_header(0); return HECMW_HECLEX_K_ABSOLUTE; }
BEAM		{ set_flag_header(0); return HECMW_HECLEX_K_BEAM; }
COMPOSITE	{ set_flag_header(0); return HECMW_HECLEX_K_COMPOSITE; }
DEFINITION	{ set_flag_header(0); return HECMW_HECLEX_K_DEFINITION; }
EGRP		{ set_flag_header(0); return HECMW_HECLEX_K_EGRP; }
GENERATE	{ set_flag_header(0); return HECMW_HECLEX_K_GENERATE; }
HECMW		{ set_flag_header(0); return HECMW_HECLEX_K_HECMW; }
INPUT		{ set_flag_header(0); return HECMW_HECLEX_K_INPUT; }
INTERFACE	{ set_flag_header(0); return HECMW_HECLEX_K_INTERFACE; }
ITEM		{ set_flag_header(0); return HECMW_HECLEX_K_ITEM; }
MATERIAL	{ set_flag_header(0); return HECMW_HECLEX_K_MATERIAL; }
MATITEM		{ set_flag_header(0); return HECMW_HECLEX_K_MATITEM; }
NAME		{ set_flag_header(0); return HECMW_HECLEX_K_NAME; }
NASTRAN		{ set_flag_header(0); return HECMW_HECLEX_K_NASTRAN; }
NGRP		{ set_flag_header(0); return HECMW_HECLEX_K_NGRP; }
NODE"-"SURF	{ set_flag_header(0); return HECMW_HECLEX_K_NODE_SURF; }
NODE"-"ELEM	{ set_flag_header(0); return HECMW_HECLEX_K_NODE_ELEM; }
RELATIVE	{ set_flag_header(0); return HECMW_HECLEX_K_RELATIVE; }
SECOPT		{ set_flag_header(0); return HECMW_HECLEX_K_SECOPT; }
SECTION		{ set_flag_header(0); return HECMW_HECLEX_K_SECTION; }
SGRP		{ set_flag_header(0); return HECMW_HECLEX_K_SGRP; }
SHELL		{ set_flag_header(0); return HECMW_HECLEX_K_SHELL; }
SOLID		{ set_flag_header(0); return HECMW_HECLEX_K_SOLID; }
STEP" "?TIME	{ set_flag_header(0); return HECMW_HECLEX_K_STEP_TIME; }
SUBITEM		{ set_flag_header(0); return HECMW_HECLEX_K_SUBITEM; }
SURF"-"SURF	{ set_flag_header(0); return HECMW_HECLEX_K_SURF_SURF; }
SYSTEM		{ set_flag_header(0); return HECMW_HECLEX_K_SYSTEM; }
TABLE		{ set_flag_header(0); return HECMW_HECLEX_K_TABLE; }
TABULAR		{ set_flag_header(0); return HECMW_HECLEX_K_TABULAR; }
TEMPERATURE	{ set_flag_header(0); return HECMW_HECLEX_K_TEMPERATURE; }
TIME		{ set_flag_header(0); return HECMW_HECLEX_K_TIME; }
TYPE		{ set_flag_header(0); return HECMW_HECLEX_K_TYPE; }
VALUE		{ set_flag_header(0); return HECMW_HECLEX_K_VALUE; }
TIMEVALUE		{ set_flag_header(0); return HECMW_HECLEX_K_TIMEVALUE; }
VALUETIME		{ set_flag_header(0); return HECMW_HECLEX_K_VALUETIME; }

=			{ set_flag_header(0); return '='; }
,			{ set_flag_header(0); return ','; }
C			{ set_flag_header(0); return 'C'; }
R			{ set_flag_header(0); return 'R'; }

{int}		{
				yylval.dval = atof(yytext);
				set_flag_header(0);
				return HECMW_HECLEX_INT;
			}
{double}	{
				yylval.dval = atof(yytext);
				set_flag_header(0);
				return HECMW_HECLEX_DOUBLE;
			}
{name}		{
				yylval.str = yytext;
				set_flag_header(0);
				return HECMW_HECLEX_NAME;
			}
{file}		{
				yylval.str = yytext;
				set_flag_header(0);
				return HECMW_HECLEX_FILENAME;
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
HECMW_heclex_get_lineno(void)
{
	return flag_including ? lineno_inc : lineno;
}


double
HECMW_heclex_get_number(void)
{
	return yylval.dval;
}


char *
HECMW_heclex_get_text(void)
{
	return yytext;
}


int
HECMW_heclex_is_including(void)
{
	return flag_including ? 1 : 0;
}


int
HECMW_heclex_next_token(void)
{
	return yylex();
}


int
HECMW_heclex_next_token_skip(int skip_token)
{
	int token;
	while((token = yylex())) {
		if(token != skip_token) break;
	}
	return token;
}


int
HECMW_heclex_set_input(FILE *fp)
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
HECMW_heclex_skip_line(void)
{
	int token;
	while((token = HECMW_heclex_next_token())) {
		if(token == HECMW_HECLEX_NL) break;
	}
	return token;
}


int
HECMW_heclex_switch_to_include(const char *filename)
{
	if(filename == NULL) return -1;
	if(flag_including) {
		HECMW_set_error(HECMW_IO_E0004, "File: %s, %s", filename);
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
HECMW_heclex_unput_token(void)
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
HECMW_hecwrap(void)
{
	return 1;
}


