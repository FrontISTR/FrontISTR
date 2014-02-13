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
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include "hecmw_config.h"
#include "hecmw_path.h"
#include "hecmw_util.h"


#ifdef __WIN32__
#define PATH_SEPARATOR '\\'
#elif __CYGWIN__
#define PATH_SEPARATOR '/'
#else
#define PATH_SEPARATOR '/'
#endif


#ifdef __WIN32__
#define isslash(c) ((c)=='/'||(c)=='\\')
#elif __CYGWIN__
#define isslash(c) ((c)=='/'||(c)=='\\')
#else
#define isslash(c) ((c)=='/')
#endif


int
HECMW_get_path_separator(void)
{
	return PATH_SEPARATOR;	
}



static int
has_drive(const char *path)
{
	const char *p = path;

	if(path == NULL) return 0;
	if(!*p || !isalpha(*p)) return 0;
	p++;
	return (*p == ':');
}


int
HECMW_is_absolute_path(const char *path)
{
	if(path == NULL) return 0;

#ifdef __WIN32__
	if(*path && isslash(*path)) return 1;
	return has_drive(path);
#elif __CYGWIN__
	if(*path && isslash(*path)) return 1;
	return has_drive(path);
#else
	return (*path && isslash(*path));
#endif
}


static char *
dir_name(const char *path)
{
	const char *p;
	static char dname[HECMW_FILENAME_LEN+1];

	if(path == NULL || strlen(path) == 0) {
		strcpy(dname, ".");
		return dname;
	}

	p = path + strlen(path) - 1;
	while(p > path && isslash(*p)) {
		p--;
	}

	while(p > path && !isslash(*p)) {
		p--;
	}

	if(p == path) {
		sprintf(dname, "%c",isslash(*p) ? PATH_SEPARATOR : '.');
		return dname;
	} else {
		do {
			p--;
		} while(p > path && isslash(*p));
	}

	if(p-path+1 > HECMW_FILENAME_LEN) {
		errno = ENAMETOOLONG;
		return NULL;
	}
	strncpy(dname, path, p-path+1);
	dname[p-path+1] = '\0';

	return dname;
}


static char *
base_name(const char *path)
{
	const char *sp,*ep;
	static char bname[HECMW_FILENAME_LEN+1];

	if(path == NULL || strlen(path) == 0) {
		strcpy(bname, ".");
		return bname;
	}

	ep = path + strlen(path) - 1;
	while(ep > path && isslash(*ep)) {
		ep--;
	}

	if(ep == path && isslash(*ep)) {
		sprintf(bname, "%c", PATH_SEPARATOR);
		return bname;
	}

	sp = ep;
	while(sp > path && !isslash(*(sp-1))) {
		sp--;
	}

	if(ep-sp+1 > HECMW_FILENAME_LEN) {
		errno = ENAMETOOLONG;
		return NULL;
	}
	strncpy(bname, sp, ep-sp+1);
	bname[ep-sp+1] = '\0';

	return bname;
}


static char *
get_bdname(const char *path, int type)
{
	const char *p, *q;
	static char bname[HECMW_FILENAME_LEN+1];
	char drive[10] = "";

	if(has_drive(path)) {
		p = path + 2;
		sprintf(drive, "%.2s", path);
	} else {
		p = path;
	}
	if(type == 'B') {	/* basename */
		q = base_name(p);
	} else {			/* dirname */
		q = dir_name(p);
	}
	if(q == NULL) return NULL;
	if(strlen(drive) > 0 && isslash(*q)) {
		if(strlen(drive) + strlen(q) > HECMW_FILENAME_LEN) {
			errno = ENAMETOOLONG;
			return NULL;
		}
		sprintf(bname, "%s%s", drive, q);
	} else {
		sprintf(bname, "%s", q);
	}
	return bname;
}


char *
HECMW_basename(const char *path)
{
	char *bname = get_bdname(path, 'B');
	if(bname == NULL) {
		HECMW_set_error(errno, "");
	}
	return bname;
}


char *
HECMW_dirname(const char *path)
{
	char *dname = get_bdname(path, 'D');
	if(dname == NULL) {
		HECMW_set_error(errno, "");
	}
	return dname;
}

