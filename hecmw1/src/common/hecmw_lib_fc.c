/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "hecmw_util.h"

char *
HECMW_strcpy_f2c(const char *fstr, int flen)
{
	int i,len;
	char *s;

	if(fstr == NULL) return NULL;
	if(flen <= 0) return NULL;

	len = 0;
	for(i=flen-1; i >= 0; i--) {
		if(fstr[i] != ' ') {
			len = i+1;
			break;
		}
	}

	if(len == 0) {
		s = HECMW_strdup("");
		if(s == NULL) {
			HECMW_set_error(errno, "");
			return NULL;
		}
		return s;
	}

	s = HECMW_malloc(len+1);
	if(s == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	strncpy(s, fstr, len);
	s[len] = '\0';
	return s;
}

char *
HECMW_strcpy_f2c_r(const char *fstr, int flen, char *buf, int bufsize)
{
	int i,len;

	if(fstr == NULL) return NULL;
	if(flen <= 0) return NULL;
	if(buf == NULL) return NULL;
	if(bufsize <= 0) return NULL; 

	len = 0;
	for(i=flen-1; i >= 0; i--) {
		if(fstr[i] != ' ') {
			len = i+1;
			break;
		}
	}
	if(len == 0) {
		buf[0] = '\0';
		return buf;
	}
	if(len > bufsize-1) {
		len = bufsize-1;
	}
	strncpy(buf, fstr, len);
	buf[len] = '\0';
	return buf;
}

int 
HECMW_strcpy_c2f(const char *cstr, char *fstr, int flen)
{
	int clen;

	if(fstr == NULL) return 0;
	if(flen <= 0) return 0;

	if(cstr == NULL) {
		clen = 0;
	} else {
		clen = strlen(cstr);
	}
	if(clen > flen) {
		clen = flen;
	}
	memset(fstr, ' ', flen);
	strncpy(fstr, cstr, clen);
	return flen;
}

