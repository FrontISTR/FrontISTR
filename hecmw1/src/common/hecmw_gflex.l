%{
#include <string.h>
#include "hecmw_gflex.h"
#include "hecmw_util.h"

typedef union {
	double dval;
	char *str;
} YYSTYPE;

static YYSTYPE yylval;
static int lineno;
static long len = 0;

%}

%option nounput 

comment_line	^(!|#).*\r?\n
comment		(!|#).*
ws			[ \t]+
nl			\r?\n

int 		[+-]?[0-9]+
double		[+-]?[0-9]*\.([0-9]*([eE][+-]?[0-9]+)?)?
name		[_a-zA-Z][_a-zA-Z0-9-]*

%%

{comment_line}	{ len += yyleng; lineno++; }
{comment}	{ len += yyleng; }
{ws}		{ len += yyleng; }
{nl}		{ len += yyleng; lineno++; return HECMW_GFLEX_NL; }

{int}		{
				len += yyleng;
				yylval.dval = atof(yytext);
				return HECMW_GFLEX_INT;
			}
{double}	{
				len += yyleng;
				yylval.dval = atof(yytext);
				return HECMW_GFLEX_DOUBLE;
			}
{name}		{
				len += yyleng;
				yylval.str = yytext;
				return HECMW_GFLEX_NAME;
			}
.			{ len += yyleng; return yytext[0]; }
<<EOF>>		{ return 0; }

%%

int
HECMW_gflex_get_lineno(void)
{
	return lineno;
}


double
HECMW_gflex_get_number(void)
{
	return yylval.dval;
}


char *
HECMW_gflex_get_text(void)
{
	return yytext;
}


int
HECMW_gflex_next_token(void)
{
	return yylex();
}


int
HECMW_gflex_next_token_skip(int skip_token)
{
	int token;
	while((token = yylex())) {
		if(token != skip_token) break;
	}
	return token;
}


long
HECMW_gflex_readsize(void)
{
	return len;
}


int
HECMW_gflex_set_input(FILE *fp)
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
	len = 0;
	return 0;
}


int
HECMW_gflex_skip_line(void)
{
	int token;
	while((token = HECMW_gflex_next_token())) {
		if(token == HECMW_GFLEX_NL) break;
	}
	return token;
}


int
HECMW_gfwrap(void)
{
	return 1;
}


