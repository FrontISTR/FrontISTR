/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
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

static int n_ptr;

static int add_info(void *ptr, size_t size, char *file, int line) {
  static struct malloc_info *info;
  int rtc;

  HECMW_assert(ptr);

  #pragma omp critical
  {
    info = malloc(sizeof(*info));
    if (info == NULL) {
      rtc = -1;
    } else {
      mem_size += size;

      info->ptr  = ptr;
      info->size = size;
      info->file = file;
      info->line = line;
      info->next = mainfo;

      mainfo = info;

      n_ptr++;
      rtc = 0;
    }
  }

  return rtc;
}

static int del_info(void *ptr) {
  struct malloc_info *p, *q;
  int rtc, i;

  HECMW_assert(ptr);

  #pragma omp critical
  {
    q = NULL;
    for (p = mainfo, i = 0; p && p->ptr != ptr; p = (q = p)->next, i++) {
      HECMW_assert(i < n_ptr);
    }
    if (p == NULL) {
      rtc = -1; /* not found */
    } else {
      if (q == NULL) {
        mainfo = p->next;
      } else {
        q->next = p->next;
      }
      mem_size -= p->size;
      free(p);

      n_ptr--;
      rtc = 0;
    }
  }
  return rtc;
}

static int change_info(void *ptrold, void *ptrnew, size_t sizenew, char *file,
                       int line) {
  struct malloc_info *p;
  long size;
  int rtc, i;

  HECMW_assert(ptrold);
  HECMW_assert(ptrnew);

  #pragma omp critical
  {
    for (p = mainfo, i = 0; p && p->ptr != ptrold; p = p->next, i++) {
      HECMW_assert(i < n_ptr);
    }
    if (p == NULL) {
      rtc = -1;
    } else {
      size = sizenew - p->size;
      mem_size += size;
      p->ptr  = ptrnew;
      p->size = sizenew;
      p->file = file;
      p->line = line;
      rtc = 0;
    }
  }
  return rtc;
}

int HECMW_list_meminfo(FILE *fp) {
  int n;
  struct malloc_info *p;

  if (fp == NULL) fp = stdout;

  n = 0;
  for (p = mainfo; p; p = p->next) {
    fprintf(fp, "HEC-MW memory info: %s:%d  ptr=%p  size=%d\n", p->file,
            p->line, p->ptr, (int)p->size);
    n++;
  }
  return n;
}

void HECMW_set_autocheck_memleak(int flag) { auto_check = flag ? 1 : 0; }

int HECMW_check_memleak(void) {
  int n;
  struct malloc_info *p;

  if (mainfo == NULL) return 0; /* no memory leaks */
  n = 0;
  for (p = mainfo; p; p = p->next) {
    fprintf(stderr,
            "HEC-MW memory check: "
            "A memory leak found at %s:%d  ptr=%p  size=%d\n",
            p->file, p->line, p->ptr, (int)p->size);
    n++;
  }
  fprintf(stderr,
          "HEC-MW memory check: "
          "%d memory leak%s found\n",
          n, (n > 1) ? "s" : "");
  return n;
}

static void check_memleak(void) { HECMW_check_memleak(); }

static int mark_check_memleak(void) {
  if (!is_check_memleak) {
    if (atexit(check_memleak) == -1) return -1;
    is_check_memleak = 1;
  }
  return 0;
}

long HECMW_get_memsize(void) { return mem_size; }

void HECMW_free_(void *ptr, char *file, int line) {
  if (ptr == NULL) return;
  if (del_info(ptr)) {
    HECMW_print_msg(HECMW_LOG_WARN, HECMW_UTIL_E9001,
                    "Not found allocated memory %p(%s:%d)\n", ptr, file, line);
  }
  free(ptr);
}

void *HECMW_malloc_(size_t size, char *file, int line) {
  void *ptr = NULL;

  ptr = malloc(size);
  if (ptr == NULL) goto error;
  if (add_info(ptr, size, file, line)) goto error;
  if (auto_check) {
    if (mark_check_memleak()) goto error;
  }
  return ptr;
error:
  free(ptr);
  return NULL;
}

void *HECMW_calloc_(size_t nmemb, size_t size, char *file, int line) {
  void *ptr = NULL;

  ptr = calloc(nmemb, size);
  if (ptr == NULL) goto error;
  if (add_info(ptr, nmemb * size, file, line)) goto error;
  if (auto_check) {
    if (mark_check_memleak()) goto error;
  }
  return ptr;
error:
  free(ptr);
  return NULL;
}

void *HECMW_realloc_(void *ptr, size_t size, char *file, int line) {
  void *ptrnew;

  ptrnew = realloc(ptr, size);

  if (size == 0 && ptr != NULL) { /* same as free */
    if (del_info(ptr)) {
      HECMW_print_msg(HECMW_LOG_WARN, HECMW_UTIL_E9001,
                      "Not found registered memory %p(%s:%d)\n", ptr, file,
                      line);
    }
    return NULL;
  }
  if (ptr == NULL) { /* same as malloc(size) */
    if (add_info(ptrnew, size, file, line)) return NULL;
  } else {
    if (ptr == ptrnew) {
      if (change_info(ptr, ptrnew, size, file, line)) {
        HECMW_print_msg(HECMW_LOG_WARN, HECMW_UTIL_E9001,
                        "Not found registered memory %p(%s:%d)\n", ptr, file,
                        line);
        if (add_info(ptrnew, size, file, line)) return NULL;
      }
    } else {
      if (del_info(ptr)) {
        HECMW_print_msg(HECMW_LOG_WARN, HECMW_UTIL_E9001,
                        "Not found registered memory %p(%s:%d)\n", ptr, file,
                        line);
      }
      if (add_info(ptrnew, size, file, line)) goto error;
    }
  }
  if (auto_check) {
    if (mark_check_memleak()) goto error;
  }
  return ptrnew;
error:
  return NULL;
}

char *HECMW_strdup_(const char *s, char *file, int line) {
  char *str = NULL;
  str       = strdup(s);
  if (str == NULL) goto error;
  if (add_info(str, strlen(str) + 1, file, line)) goto error;
  if (auto_check) {
    if (mark_check_memleak()) goto error;
  }
  return str;
error:
  free(str);
  return NULL;
}

#define SMALL_BUFFER_SIZE 240
#define BUFFER_NUM 256
#define POOL_SIZE 64
void *small_buffer_pool[POOL_SIZE];
int small_buffer_pool_cap = 0;

void *HECMW_malloc_gpu(size_t size) {
  if(size <= SMALL_BUFFER_SIZE) {
    char *r;
    if(small_buffer_pool_cap > 0) {
      r = small_buffer_pool[small_buffer_pool_cap - 1];
      small_buffer_pool_cap--;
    } else {
      static char *shared_buffer = NULL;
      static int shared_buffer_used;
      if(shared_buffer == NULL) {
        shared_buffer = (char*)malloc((SMALL_BUFFER_SIZE + 16) * BUFFER_NUM);
        shared_buffer_used = 0;
        for(int i = 0; i < BUFFER_NUM; i++) {
          *(size_t*)(shared_buffer + (SMALL_BUFFER_SIZE + 16) * i) = i + 1;
        }
      }

      r = shared_buffer + (SMALL_BUFFER_SIZE + 16) * shared_buffer_used;

      shared_buffer_used++;
      if(shared_buffer_used >= BUFFER_NUM) {
        shared_buffer = NULL;
      }
    }

    *(size_t*)(r + 8) = size;
    return r + 16;
  } else {
    void *p = malloc(size + 8);
    *(size_t*)p = size;
    return (void*)((char*)p + 8);
  }
}

void HECMW_free_gpu(void *ptr) {
  if(ptr != NULL) {
    int old_size = *(size_t*)((char*)ptr - 8);
    if(old_size <= SMALL_BUFFER_SIZE) {
      if(small_buffer_pool_cap >= POOL_SIZE) {
        void *free_ptr = small_buffer_pool[0];
        for(int i = 0; i < POOL_SIZE - 1; i++) {
          small_buffer_pool[i] = small_buffer_pool[i + 1];
        }
        small_buffer_pool_cap--;

        int offset = *(size_t*)free_ptr - 1;
        *(size_t*)free_ptr = 0;
        char *base = (char*)free_ptr - (SMALL_BUFFER_SIZE + 16) * offset;

        bool in_use = false;
        for(int i = 0; i < BUFFER_NUM; i++) {
          if(*(size_t*)(base + (SMALL_BUFFER_SIZE + 16) * i) != 0) {
            in_use = true;
            break;
          }
        }

        if(!in_use) {
          free(base);
        }
      }

      small_buffer_pool[small_buffer_pool_cap] = (char*)ptr - 16;
      small_buffer_pool_cap++;
    } else {
      free((void*)((char*)ptr - 8));
    }
  }
}

void *HECMW_calloc_gpu(size_t nmemb, size_t size) {
  void *p = HECMW_malloc_gpu(nmemb * size);
  memset(p, 0, nmemb * size);
  return p;
}

void *HECMW_realloc_gpu(void *ptr, size_t size) {
  void *p = HECMW_malloc_gpu(size);
  if(ptr != NULL) {
    int old_size = *(size_t*)((char*)ptr - 8);
    int copy_size = old_size;
    if(copy_size >= size) {
      copy_size = size;
    }
    memcpy(p, ptr, copy_size);
    HECMW_free_gpu(ptr);
  }
  return p;
}

char *HECMW_strdup_gpu(const char *s) {
  int l = strlen(s);
  void *p = HECMW_malloc_gpu(l + 1);
  memcpy(p, (void*)s, l + 1);
  return (char*)p;
}
