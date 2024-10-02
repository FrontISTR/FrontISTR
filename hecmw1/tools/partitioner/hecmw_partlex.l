%{
#include <string.h>
#include <errno.h>
#include "hecmw_partlex.h"

typedef union {
  double dval;
  char *str;
} YYSTYPE;

static YYSTYPE yylval;
static int lineno;
static long len = 0;

%}

%option nounput

comment		^(!!|#).*\r?\n
ws		[ \t]+
nl		\r?\n

int 		[+-]?[0-9]+
double		[+-]?[0-9]*\.([0-9]*([eE][+-]?[0-9]+)?)?
name		[_a-zA-Z][_a-zA-Z0-9-]*
file		[/._a-zA-Z0-9-]+

%%

{comment}	{ len += yyleng; lineno++; }
{ws}		{ len += yyleng; }
{nl}		{ len += yyleng; lineno++; return HECMW_PARTLEX_NL; }

^!PARTITION	{ len += yyleng; return HECMW_PARTLEX_H_PARTITION; }

METHOD		{ len += yyleng; return HECMW_PARTLEX_K_METHOD; }
TYPE		{ len += yyleng; return HECMW_PARTLEX_K_TYPE; }
DOMAIN		{ len += yyleng; return HECMW_PARTLEX_K_DOMAIN; }
DEPTH		{ len += yyleng; return HECMW_PARTLEX_K_DEPTH; }
UCD		{ len += yyleng; return HECMW_PARTLEX_K_UCD; }
CONTACT		{ len += yyleng; return HECMW_PARTLEX_K_CONTACT; }
PART		{ len += yyleng; return HECMW_PARTLEX_K_PART; }

RCB		{ len += yyleng; return HECMW_PARTLEX_V_RCB; }
KMETIS		{ len += yyleng; return HECMW_PARTLEX_V_KMETIS; }
PMETIS		{ len += yyleng; return HECMW_PARTLEX_V_PMETIS; }
NODE-BASED	{ len += yyleng; return HECMW_PARTLEX_V_NODE_BASED; }
USER		{ len += yyleng; return HECMW_PARTLEX_V_USER; }
ELEMENT-BASED	{ len += yyleng; return HECMW_PARTLEX_V_ELEMENT_BASED; }
DEFAULT		{ len += yyleng; return HECMW_PARTLEX_V_DEFAULT; }
AGGREGATE	{ len += yyleng; return HECMW_PARTLEX_V_AGGREGATE; }
DISTRIBUTE	{ len += yyleng; return HECMW_PARTLEX_V_DISTRIBUTE; }
SIMPLE		{ len += yyleng; return HECMW_PARTLEX_V_SIMPLE; }

=		{ len += yyleng; return '='; }
,		{ len += yyleng; return ','; }
x		{ len += yyleng; return 'x'; }
y		{ len += yyleng; return 'y'; }
z		{ len += yyleng; return 'z'; }

{int}		{
			len += yyleng;
			yylval.dval = atof(yytext);
			return HECMW_PARTLEX_INT;
		}
{double}	{
			len += yyleng;
			yylval.dval = atof(yytext);
			return HECMW_PARTLEX_DOUBLE;
		}
{name}		{
			len += yyleng;
			yylval.str = yytext;
			return HECMW_PARTLEX_NAME;
		}
{file}		{
			len += yyleng;
			yylval.str = yytext;
			return HECMW_PARTLEX_FILENAME;
		}
.		{ len += yyleng; return yytext[0]; }
<<EOF>>		{ return 0; }

%%

int
HECMW_partlex_get_lineno(void)
{
  return lineno;
}


double
HECMW_partlex_get_number(void)
{
  return yylval.dval;
}


char *
HECMW_partlex_get_text(void)
{
  return yytext;
}


int
HECMW_partlex_next_token(void)
{
  return yylex();
}


int
HECMW_partlex_next_token_skip(int skip_token)
{
  int token;
  while((token = yylex())) {
    if(token != skip_token) break;
  }
  return token;
}


long
HECMW_partlex_readsize(void)
{
  return len;
}


int
HECMW_partlex_set_input(FILE *fp)
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
HECMW_partlex_skip_line(void)
{
  int token;
  while((token = HECMW_partlex_next_token())) {
    if(token == HECMW_PARTLEX_NL) break;
  }
  return token;
}


int
HECMW_partlex_unput_token(void)
{
  char *p = yytext;
  while(*p) {
    if(*p == '\n') lineno--;
    len--;
    p++;
  }
  yyless(0);
  return 0;
}


int
HECMW_partwrap(void)
{
  return 1;
}
