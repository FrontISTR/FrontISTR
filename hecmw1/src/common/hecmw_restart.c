/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_config.h"
#include "hecmw_restart.h"

struct restart_list {
  void *data;
  size_t size;
  struct restart_list *next;
};

static struct restart_list *first_list = NULL;
static struct restart_list *previ_list = NULL;
static FILE *restart_fp;

static void clear(void) {
  struct restart_list *p, *q;

  for (p = first_list; p; p = q) {
    q = p->next;
    HECMW_free(p->data);
    HECMW_free(p);
  }
  first_list = NULL;
  previ_list = NULL;
}

int HECMW_restart_open_by_name(char *name_ID) {
  char *filename;

  if (name_ID) {
    if ((filename = HECMW_ctrl_get_restart_file(name_ID)) == NULL) return -1;
  } else {
    /* io is bitmap */
    int io = HECMW_CTRL_FILE_IO_IN | HECMW_CTRL_FILE_IO_INOUT;
    if ((filename = HECMW_ctrl_get_restart_file_by_io(io)) == NULL) return -1;
  }

  if ((restart_fp = fopen(filename, "rb")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0101, "File: %s, %s", filename,
                    HECMW_strmsg(errno));
    HECMW_free(filename);
    return -1;
  }
  HECMW_free(filename);
  return 0;
}

int HECMW_restart_open(void) { return HECMW_restart_open_by_name(NULL); }

int HECMW_restart_close(void) {
  if (restart_fp == NULL) return 0;
  if (fclose(restart_fp)) {
    HECMW_set_error(HECMW_UTIL_E0102, HECMW_strmsg(errno));
    return -1;
  }
  restart_fp = NULL;
  return 0;
}

static int restart_read(void *ptr, size_t size, size_t nmemb, FILE *stream) {
  int rc, msgno;

  rc = fread(ptr, size, nmemb, stream);
  if (rc != nmemb) {
    if (feof(stream)) {
      msgno = HECMW_UTIL_E0104;
    } else if (ferror(stream)) {
      msgno = HECMW_UTIL_E0105;
    }
    return msgno;
  }
  return 0;
}

void *HECMW_restart_read(void *addr) {
  int rc;
  size_t size;
  void *data;

  if (restart_fp == NULL) {
    HECMW_set_error(HECMW_UTIL_E0103, "");
    return NULL;
  }

  /* read size */
  rc = restart_read(&size, sizeof(size), 1, restart_fp);
  if (rc) {
    HECMW_set_error(rc, "");
    return NULL;
  }

  /* read data */
  if (addr == NULL) {
    data = HECMW_malloc(size);
    if (data == NULL) {
      HECMW_set_error(errno, "");
      return NULL;
    }
  } else {
    data = addr;
  }
  rc = restart_read(data, size, 1, restart_fp);
  if (rc) {
    HECMW_set_error(rc, "");
    return NULL;
  }

  return data;
}

int HECMW_restart_add(void *data, size_t size) {
  struct restart_list *rst;

  rst = HECMW_malloc(sizeof(*rst));
  if (rst == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  rst->data = data;
  rst->size = size;
  rst->next = NULL;

  if (previ_list != NULL) previ_list->next = rst;
  previ_list                               = rst;
  if (first_list == NULL) first_list       = rst;

  return 0;
}

int HECMW_restart_add_int(int *data, int n_data) {
  return HECMW_restart_add(data, sizeof(int) * n_data);
}

int HECMW_restart_add_double(double *data, int n_data) {
  return HECMW_restart_add(data, sizeof(double) * n_data);
}

static int restart_write(const void *ptr, size_t size, size_t nmemb,
                         FILE *stream) {
  int rc, msgno;

  rc = fwrite(ptr, size, nmemb, stream);
  if (rc != nmemb) {
    if (feof(stream)) {
      msgno = HECMW_UTIL_E0104;
    } else if (ferror(stream)) {
      msgno = HECMW_UTIL_E0105;
    }
    return msgno;
  }
  return 0;
}

int HECMW_restart_write_by_name(char *name_ID) {
  FILE *fp;
  struct restart_list *p;
  char *filename;
  int rc;

  if (previ_list == NULL) return 0;

  if (name_ID) {
    if ((filename = HECMW_ctrl_get_restart_file(name_ID)) == NULL) return -1;
  } else {
    /* io is bitmap */
    int io = HECMW_CTRL_FILE_IO_OUT | HECMW_CTRL_FILE_IO_INOUT;
    if ((filename = HECMW_ctrl_get_restart_file_by_io(io)) == NULL) return -1;
  }

  if (HECMW_ctrl_is_subdir()) {
    if (HECMW_ctrl_make_subdir(filename)) {
      HECMW_free(filename);
      return -1;
    }
  }

  if ((fp = fopen(filename, "wb")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0101, "File: %s, %s", filename,
                    HECMW_strmsg(errno));
    HECMW_free(filename);
    return -1;
  }
  HECMW_free(filename);

  for (p = first_list; p; p = p->next) {
    /* write size */
    rc = restart_write(&p->size, sizeof(p->size), 1, fp);
    if (rc) {
      HECMW_set_error(rc, "");
      return -1;
    }
    /* write data */
    rc = restart_write(p->data, p->size, 1, fp);
    if (rc) {
      HECMW_set_error(rc, "");
      return -1;
    }
  }

  if (fclose(fp)) {
    HECMW_set_error(errno, "");
    return -1;
  }

  clear();

  return 0;
}

int HECMW_restart_write(void) { return HECMW_restart_write_by_name(NULL); }

/*----------------------------------------------------------------------------*/

void hecmw_restart_open_by_name_if(char *name_ID, int *err, int len) {
  char *name = NULL;
  char cname[HECMW_NAME_LEN + 1];

  *err = 1;

  if (name_ID) {
    if (HECMW_strcpy_f2c_r(name_ID, len, cname, sizeof(cname)) == NULL) return;
    name = cname;
  }

  if (HECMW_restart_open_by_name(name)) return;

  *err = 0;
}

void hecmw_restart_open_by_name_if_(char *name_ID, int *err, int len) {
  hecmw_restart_open_by_name_if(name_ID, err, len);
}

void hecmw_restart_open_by_name_if__(char *name_ID, int *err, int len) {
  hecmw_restart_open_by_name_if(name_ID, err, len);
}

void HECMW_RESTART_OPEN_BY_NAME_IF(char *name_ID, int *err, int len) {
  hecmw_restart_open_by_name_if(name_ID, err, len);
}

/*----------------------------------------------------------------------------*/

void hecmw_restart_open_if(int *err) {
  hecmw_restart_open_by_name_if(NULL, err, 0);
}

void hecmw_restart_open_if_(int *err) { hecmw_restart_open_if(err); }

void hecmw_restart_open_if__(int *err) { hecmw_restart_open_if(err); }

void HECMW_RESTART_OPEN_IF(int *err) { hecmw_restart_open_if(err); }

/*----------------------------------------------------------------------------*/

void hecmw_restart_close_if(int *err) { *err = HECMW_restart_close() ? 1 : 0; }

void hecmw_restart_close_if_(int *err) { hecmw_restart_close_if(err); }

void hecmw_restart_close_if__(int *err) { hecmw_restart_close_if(err); }

void HECMW_RESTART_CLOSE_IF(int *err) { hecmw_restart_close_if(err); }

/*----------------------------------------------------------------------------*/

void hecmw_restart_read_int_if(int *dst, int *err) {
  *err = 1;

  if (dst == NULL) {
    HECMW_set_error(HECMW_ALL_E0101, "");
    return;
  }
  if (HECMW_restart_read(dst) == NULL) return;

  *err = 0;
}

void hecmw_restart_read_int_if_(int *dst, int *err) {
  hecmw_restart_read_int_if(dst, err);
}

void hecmw_restart_read_int_if__(int *dst, int *err) {
  hecmw_restart_read_int_if(dst, err);
}

void HECMW_RESTART_READ_INT_IF(int *dst, int *err) {
  hecmw_restart_read_int_if(dst, err);
}

/*----------------------------------------------------------------------------*/

void hecmw_restart_read_real_if(double *dst, int *err) {
  *err = 1;

  if (dst == NULL) {
    HECMW_set_error(HECMW_ALL_E0101, "");
    return;
  }
  if (HECMW_restart_read(dst) == NULL) return;

  *err = 0;
}

void hecmw_restart_read_real_if_(double *dst, int *err) {
  hecmw_restart_read_real_if(dst, err);
}

void hecmw_restart_read_real_if__(double *dst, int *err) {
  hecmw_restart_read_real_if(dst, err);
}

void HECMW_RESTART_READ_REAL_IF(double *dst, int *err) {
  hecmw_restart_read_real_if(dst, err);
}

/*----------------------------------------------------------------------------*/

static void *restart_add_alloc(void *data, int byte, int n_data) {
  int size;
  void *cdata;

  size  = byte * n_data;
  cdata = HECMW_malloc(size);
  if (cdata == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }
  memcpy(cdata, data, size);

  return cdata;
}

void hecmw_restart_add_int_if(int *data, int *n_data, int *err) {
  int *cdata;

  cdata = restart_add_alloc(data, sizeof(int), *n_data);
  if (cdata == NULL) {
    *err = 1;
    return;
  }

  *err = HECMW_restart_add_int(cdata, *n_data);
}

void hecmw_restart_add_int_if_(int *data, int *n_data, int *err) {
  hecmw_restart_add_int_if(data, n_data, err);
}

void hecmw_restart_add_int_if__(int *data, int *n_data, int *err) {
  hecmw_restart_add_int_if(data, n_data, err);
}

void HECMW_RESTART_ADD_INT_IF(int *data, int *n_data, int *err) {
  hecmw_restart_add_int_if(data, n_data, err);
}

/*----------------------------------------------------------------------------*/

void hecmw_restart_add_real_if(double *data, int *n_data, int *err) {
  void *cdata;

  cdata = restart_add_alloc(data, sizeof(double), *n_data);
  if (cdata == NULL) {
    *err = 1;
    return;
  }
  *err = HECMW_restart_add_double(cdata, *n_data);
}

void hecmw_restart_add_real_if_(double *data, int *n_data, int *err) {
  hecmw_restart_add_real_if(data, n_data, err);
}

void hecmw_restart_add_real_if__(double *data, int *n_data, int *err) {
  hecmw_restart_add_real_if(data, n_data, err);
}

void HECMW_RESTART_ADD_REAL_IF(double *data, int *n_data, int *err) {
  hecmw_restart_add_real_if(data, n_data, err);
}

/*----------------------------------------------------------------------------*/

void hecmw_restart_write_by_name_if(char *name_ID, int *err, int len) {
  char *name = NULL;
  char cname[HECMW_NAME_LEN + 1];

  *err = 1;

  if (name_ID) {
    if (HECMW_strcpy_f2c_r(name_ID, len, cname, sizeof(cname)) == NULL) return;
    name = cname;
  }

  if (HECMW_restart_write_by_name(name)) return;

  *err = 0;
}

void hecmw_restart_write_by_name_if_(char *name_ID, int *err, int len) {
  hecmw_restart_write_by_name_if(name_ID, err, len);
}

void hecmw_restart_write_by_name_if__(char *name_ID, int *err, int len) {
  hecmw_restart_write_by_name_if(name_ID, err, len);
}

void HECMW_RESTART_WRITE_BY_NAME_IF(char *name_ID, int *err, int len) {
  hecmw_restart_write_by_name_if(name_ID, err, len);
}

/*----------------------------------------------------------------------------*/

void hecmw_restart_write_if(int *err) {
  hecmw_restart_write_by_name_if(NULL, err, 0);
}

void hecmw_restart_write_if_(int *err) { hecmw_restart_write_if(err); }

void hecmw_restart_write_if__(int *err) { hecmw_restart_write_if(err); }

void HECMW_RESTART_WRITE_IF(int *err) { hecmw_restart_write_if(err); }
