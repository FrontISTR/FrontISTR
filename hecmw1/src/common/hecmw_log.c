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
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include "hecmw_config.h"
#include "hecmw_util.h"
#include "hecmw_log.h"


#define HECMW_LINE_MAX 1023


#define USEFUL_LOGLVS(lv) ((lv) ? (lv|(lv-1)) : 0)


static struct loglv_ent {
	int loglv;	
	char *str;	
} hecmw_loglv_table[] = {
	{HECMW_LOG_ERROR,"Error"},{HECMW_LOG_WARN, "Warning"},
	{HECMW_LOG_INFO, "Info" },{HECMW_LOG_DEBUG,"Debug"},
};


struct log_ent {
	FILE *fp;	
	char file[HECMW_FILENAME_LEN+1];	
	int lv;		
	int opt;	
};


static struct log_ent logent[HECMW_LOG_MAX];


static int useful_logent[HECMW_LOG_MAX];


static int loglevels = USEFUL_LOGLVS(HECMW_LOG_INFO);


static int enable_by_rank = 1;


static int
check_loglv(int loglv)
{
	int i;
	for(i=0; i < sizeof(hecmw_loglv_table)/sizeof(hecmw_loglv_table[0]); i++) {
		if(loglv == hecmw_loglv_table[i].loglv) return 0;
	}
	return -1;	/* not found */
}

void
HECMW_setloglv(int loglv)
{
	if(loglv != HECMW_LOG_NONE && check_loglv(loglv)) {
		return;
	}
	loglevels = USEFUL_LOGLVS(loglv);
}

int
HECMW_openlog(const char *logfile, int loglv, int options)
{
	FILE *fp;
	int i,newent_idx;
	struct log_ent *newent;
	char rank[10];
	char logfilename[HECMW_FILENAME_LEN+1];

	if(logfile == NULL) {
		HECMW_set_error(HECMW_UTIL_E9011, "Not specified log filename");
		return -1;
	}
	HECMW_snprintf(rank, sizeof(rank), ".%d", HECMW_comm_get_rank());
	if((strlen(logfile) + strlen(rank)) > HECMW_FILENAME_LEN) {
		HECMW_set_error(HECMW_UTIL_E9011, "Filename too long");
		return -1;
	}
	sprintf(logfilename, "%s%s", logfile, rank);

	loglv &= HECMW_LOG_ALL;
	if(!loglv) {
		HECMW_set_error(HECMW_UTIL_E9011, "Invalid log level");
		return -1;
	}

	options &= HECMW_LOG_OPTALL;

	/* search same logfile */
	for(i=0; i < HECMW_LOG_MAX; i++) {
		if(!useful_logent[i]) continue;
		if(strcmp(logfilename, logent[i].file) == 0) {
			logent[i].lv = loglv;		/* override */
			logent[i].opt = options;	/* override */
			return 0;
		}
	}
	/* search free entry */
	newent = NULL;
	for(i=0; i < HECMW_LOG_MAX; i++) {
		if(useful_logent[i]) continue;	/* using */
		newent = &logent[i];
		newent_idx = i;
		break;
	}
	if(newent == NULL) {	/* no free entry */
		HECMW_set_error(HECMW_UTIL_E9011, "No free entry");
		return -1;
	}
	/* regist new log entry */
	strcpy(newent->file, logfilename);
	newent->lv = loglv;
	newent->opt = options;
	if(enable_by_rank && loglv & loglevels) {
		if((fp = fopen(logfilename, "a")) == NULL) {
			HECMW_set_error(HECMW_UTIL_E9011, "File %s, %s", logfilename, strerror(errno));
			return -1;
		}
	} else {
		fp = NULL;
	}
	newent->fp = fp;
	/* mark */
	useful_logent[newent_idx] = 1;
	return newent_idx + 1;
}


int
HECMW_closelog(int id)
{
	struct log_ent *p; 

	if(id <= 0 || id > HECMW_LOG_MAX) {
		HECMW_set_error(HECMW_UTIL_E9013, "No such log file");
		return -1;
	}

	p = (struct log_ent*)&logent[id];
	if(p->fp == NULL) {
		if(fclose(p->fp)) {
			HECMW_set_error(HECMW_UTIL_E9013, "File %s, %s", p->file, strerror(errno));
			return -1;
		}
	}
	memset(p, 0, sizeof(*p));
	useful_logent[id] = 0;

	return 0;
}


static void
output_log(int loglv, const char *fmt, va_list ap, FILE *fp)
{
	int i,len;
	char *p;
	char buf[HECMW_LINE_MAX+1];

	HECMW_assert(fp);

	/* date */
	p = HECMW_get_date();
	if(p == NULL) {
		p = "Could not get date";
	}
	strcpy(buf, p);

	/* log level */
	for(i=0; i < sizeof(hecmw_loglv_table)/sizeof(hecmw_loglv_table[0]); i++) {
		if(hecmw_loglv_table[i].loglv == loglv) {
			p = hecmw_loglv_table[i].str;
			break;
		}
	}
	len = strlen(buf);
	HECMW_snprintf(buf+len, sizeof(buf)-len, " %s: ", p);

	/* output header */
	fprintf(fp, "%s", buf);

	if(fmt == NULL) {
		fprintf(fp, "\n");
		return;
	}

	HECMW_vsnprintf(buf, sizeof(buf), fmt, ap);		/* message body */

	/* output body */
	fputs(buf, fp);

	/* add '\n' if needed */
	len = strlen(buf);
	if((len > 0 && buf[len-1] != '\n') || len == 0) {
		fprintf(fp, "\n");
	}

	fflush(fp);
}


int
HECMW_vlog(int loglv, const char *fmt, va_list ap)
{
	int i, is_output;

	if (!enable_by_rank) return 0;
	if(check_loglv(loglv)) {
		HECMW_set_error(HECMW_UTIL_E9012, "Invalid log level");
		return -1;
	}

	if(!(loglv & loglevels)) return 0;

	is_output = 0;
	for(i=0; i < HECMW_LOG_MAX; i++) {
		struct log_ent *p = &logent[i];
		if(!useful_logent[i]) continue;
		if(loglv & loglevels & p->lv) {
			if(p->fp == NULL) {
				if((p->fp = fopen(p->file, "a")) == NULL) {
					HECMW_set_error(HECMW_UTIL_E9011, "File %s, %s", p->file, strerror(errno));
					return -1;
				}
			}
			output_log(loglv, fmt, ap, p->fp);
			if(p->opt & HECMW_LOG_PERROR) {
				output_log(loglv, fmt, ap, stderr);
			}
			is_output = 1;
		}
	}
	if(!is_output) {
		output_log(loglv, fmt, ap, stderr);
	}

	return 0;
}


int
HECMW_log(int loglv, const char *fmt, ...)
{
	int rc;
	va_list ap;
	va_start(ap, fmt);
	rc = HECMW_vlog(loglv, fmt, ap);
	va_end(ap);
	return rc;
}


void
HECMW_log_set_enable(int from, int to, int true_or_false)
{
	if (from > to) return;
	if (from <= HECMW_comm_get_rank() && HECMW_comm_get_rank() <= to) {
		enable_by_rank = (true_or_false) ? 1 : 0;
	}
}


 


void
hecmw_log_if(int *loglv, char *msg, int len)
{
	char s[HECMW_MSG_LEN+1];

	if(len > HECMW_MSG_LEN) {
		len = HECMW_MSG_LEN;
	}

	if(HECMW_strcpy_f2c_r(msg, len, s, sizeof(s)) == NULL) return;
	HECMW_log(*loglv, s);
}



void
hecmw_log_if_(int *loglv, char *msg, int len)
{
	hecmw_log_if(loglv, msg, len);
}



void
hecmw_log_if__(int *loglv, char *msg, int len)
{
	hecmw_log_if(loglv, msg, len);
}



void
HECMW_LOG_IF(int *loglv, char *msg, int len)
{
	hecmw_log_if(loglv, msg, len);
}


/*----------------------------------------------------------------------------*/


void
hecmw_setloglv_if(int *loglv)
{
	HECMW_setloglv(*loglv);
}



void
hecmw_setloglv_if_(int *loglv)
{
	hecmw_setloglv_if(loglv);
}



void
hecmw_setloglv_if__(int *loglv)
{
	hecmw_setloglv_if(loglv);
}



void
HECMW_SETLOGLV_IF(int *loglv)
{
	hecmw_setloglv_if(loglv);
}


/*----------------------------------------------------------------------------*/


void
hecmw_log_set_enable_if(int *from, int *to, int *true_or_false)
{
	HECMW_log_set_enable(*from, *to, *true_or_false);
}



void
hecmw_log_set_enable_if_(int *from, int *to, int *true_or_false)
{
	hecmw_log_set_enable_if(from, to, true_or_false);
}



void
hecmw_log_set_enable_if__(int *from, int *to, int *true_or_false)
{
	hecmw_log_set_enable_if(from, to, true_or_false);
}



void
HECMW_LOG_SET_ENABLE_IF(int *from, int *to, int *true_or_false)
{
	hecmw_log_set_enable_if(from, to, true_or_false);
}


/*----------------------------------------------------------------------------*/


void
hecmw_openlog_if(char *logfile, int *loglv, int *options, int *id, int *ierror, int len)
{
	char buf[HECMW_NAME_LEN+1];
	if(HECMW_strcpy_f2c_r(logfile, len, buf, sizeof(buf)) == NULL) {
		*ierror = 1;
		return;
	}
	if((*id = HECMW_openlog(buf, *loglv, *options)) == -1) {
		*ierror = 1;
		return;
	}
	*ierror = 0;
}



void
hecmw_openlog_if_(char *logfile, int *loglv, int *options, int *id, int *ierror, int len)
{
	hecmw_openlog_if(logfile, loglv, options, id, ierror, len);
}



void
hecmw_openlog_if__(char *logfile, int *loglv, int *options, int *id, int *ierror, int len)
{
	hecmw_openlog_if(logfile, loglv, options, id, ierror, len);
}



void
HECMW_OPENLOG_IF(char *logfile, int *loglv, int *options, int *id, int *ierror, int len)
{
	hecmw_openlog_if(logfile, loglv, options, id, ierror, len);
}


/*----------------------------------------------------------------------------*/


void
hecmw_closelog_if(int *id, int *ierror)
{
	if(HECMW_closelog(*id)) {
		*ierror = 1;
		return;
	}
	*ierror = 0;
}



void
hecmw_closelog_if_(int *id, int *ierror)
{
	hecmw_closelog_if(id, ierror);
}



void
hecmw_closelog_if__(int *id, int *ierror)
{
	hecmw_closelog_if(id, ierror);
}



void
HECMW_CLOSELOG_IF(int *id, int *ierror)
{
	hecmw_closelog_if(id, ierror);
}
