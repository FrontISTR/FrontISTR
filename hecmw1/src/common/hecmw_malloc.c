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



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "hecmw_config.h"
#include "hecmw_util.h"
#include "hecmw_malloc.h"


struct malloc_info {
	void *ptr;		
	size_t size;	
	char *file;		
	int line;		
	struct malloc_info *next;	
};


static struct malloc_info *mainfo;


static int is_check_memleak;


static long mem_size;


static int auto_check = 1;



static int 
add_info(void *ptr, size_t size, char *file, int line)
{
	static struct malloc_info *info;

	HECMW_assert(ptr);

	info = malloc(sizeof(*info));
	if(info == NULL) return -1;

	mem_size += size;

	info->ptr = ptr;
	info->size = size;
	info->file = file;
	info->line = line;
	info->next = mainfo;

	mainfo = info;

	return 0;
}



static int 
del_info(void *ptr)
{
	struct malloc_info *p,*q;

	HECMW_assert(ptr);

	q = NULL;
	for(p=mainfo; p && p->ptr!=ptr; p=(q=p)->next);
	if(p == NULL) return -1;  /* not found */

	if(q == NULL) {
		mainfo = p->next;
	} else {
		q->next = p->next;
	}
	mem_size -= p->size;
	free(p);
	return 0;
}



static int 
change_info(void *ptrold, void *ptrnew, size_t sizenew, char *file, int line)
{
	struct malloc_info *p;
	size_t size;

	HECMW_assert(ptrold);
	HECMW_assert(ptrnew);

	for(p=mainfo; p && p->ptr!=ptrold; p=p->next);
	if(p == NULL) return -1;

	size = sizenew - p->size;
	mem_size += size;
	p->ptr = ptrnew;
	p->size = sizenew;
	p->file = file;
	p->line = line;
	return 0;
}


int
HECMW_list_meminfo(FILE *fp)
{
	int n;
	struct malloc_info *p;

	if(fp == NULL) fp = stdout; 

	n = 0;
	for(p=mainfo; p; p=p->next) {
		fprintf(fp, "HEC-MW memory info: %s:%d  ptr=%p  size=%d\n",
					p->file, p->line, p->ptr, (int)p->size);
		n++;
	}
	return n;
}


void
HECMW_set_autocheck_memleak(int flag)
{
	auto_check = flag ? 1 : 0;
}


int
HECMW_check_memleak(void)
{
	int n;
	struct malloc_info *p;

	if(mainfo == NULL) return 0;    /* no memory leaks */
	n = 0;
	for(p=mainfo; p; p=p->next) {
		fprintf(stderr, "HEC-MW memory check: "
				"A memory leak found at %s:%d  ptr=%p  size=%d\n",
				p->file, p->line, p->ptr, (int)p->size);
		n++;
	}
	fprintf(stderr, "HEC-MW memory check: "
			"%d memory leak%s found\n", n, (n > 1) ? "s" : "");
	return n;
}



static void
check_memleak(void)
{
	HECMW_check_memleak();
}



static int
mark_check_memleak(void)
{
	if(!is_check_memleak) {
		if(atexit(check_memleak) == -1) return -1;
		is_check_memleak = 1;
	}
	return 0;
}


long
HECMW_get_memsize(void)
{
	return mem_size;
}


void
HECMW_free_(void *ptr, char *file, int line)
{
	if(ptr == NULL) return;
	if(del_info(ptr)) {
		HECMW_print_msg(HECMW_LOG_WARN, HECMW_UTIL_E9001,
			"Not found allocated memory %p(%s:%d)\n", ptr,file,line);
	}
	free(ptr);
}


void *
HECMW_malloc_(size_t size, char *file, int line) 
{
	void *ptr = NULL;
	
	ptr = malloc(size);
	if(ptr == NULL) goto error;
	if(add_info(ptr, size, file, line)) goto error;
	if(auto_check) {
		if(mark_check_memleak()) goto error;
	}
	return ptr;
error:
	free(ptr);
	return NULL;
}


void *
HECMW_calloc_(size_t nmemb, size_t size, char *file, int line) 
{
	void *ptr = NULL;
	
	ptr = calloc(nmemb, size);
	if(ptr == NULL) goto error;
	if(add_info(ptr, nmemb * size, file, line)) goto error;
	if(auto_check) {
		if(mark_check_memleak()) goto error;
	}
	return ptr;
error:
	free(ptr);
	return NULL;
}


void *
HECMW_realloc_(void *ptr, size_t size, char *file, int line) 
{
	void *ptrnew;

	ptrnew = realloc(ptr, size);

	if(size == 0 && ptr != NULL) {	/* same as free */
		if(del_info(ptr)) {
			HECMW_print_msg(HECMW_LOG_WARN, HECMW_UTIL_E9001,
				"Not found registered memory %p(%s:%d)\n",
				ptr,file,line);
		}
		return NULL;
	}
	if(ptr == NULL) {    /* same as malloc(size) */
		if(add_info(ptrnew, size, file, line)) return NULL;
	} else {
		if(ptr == ptrnew) {
			if(change_info(ptr, ptrnew, size, file, line)) {
				HECMW_print_msg(HECMW_LOG_WARN, HECMW_UTIL_E9001,
					"Not found registered memory %p(%s:%d)\n",
					ptr,file,line);
				if(add_info(ptrnew, size, file, line)) return NULL;
			}
		} else {
			if(del_info(ptr)) {
				HECMW_print_msg(HECMW_LOG_WARN, HECMW_UTIL_E9001,
						"Not found registered memory %p(%s:%d)\n", ptr,file,line);
			}
			if(add_info(ptrnew, size, file, line)) goto error;
		}
	}
	if(auto_check) {
		if(mark_check_memleak()) goto error;
	}
	return ptrnew;
error:
	return NULL;
}


char *
HECMW_strdup_(const char *s, char *file, int line)
{
	char *str = NULL;
	str = strdup(s);
	if(str == NULL) goto error;
	if(add_info(str, strlen(str)+1, file, line)) goto error;
	if(auto_check) {
		if(mark_check_memleak()) goto error;
	}
	return str;
error:
	free(str);
	return NULL;
}
