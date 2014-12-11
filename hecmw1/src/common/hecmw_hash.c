/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *                                                                     *
 *     Last Update : 2014/07/29                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Naoki MORITA (GSFS, the Univ. of Tokyo)       *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/
 
#include "hecmw_hash.h"
#include "hecmw_io_struct.h"

typedef struct List List;
typedef struct Bin Bin;

int hecmw_hash_p_resize(hecmw_hash_p *hash);
static unsigned int hash_key(const char *str);
static List *get_list(Bin *bin, const char *key);
static unsigned int hecmw_hashsize_init_index = 0;

struct List {
  unsigned int hashkey;
  char *key;
  void *value;
};

struct Bin {
  unsigned int n;
  List *list;
};

struct hecmw_hash_p {
  unsigned int n;
  unsigned int num_put;
  Bin *bin;
};

hecmw_hash_p *hash_ng; /* node group */
hecmw_hash_p *hash_eg; /* element group */
hecmw_hash_p *hash_sg; /* surface group */

static const unsigned int hecmw_hashsize_template[] =
{
  1021,
  2039,
  4093,
  8191,
  16381,
  32749,
  65521,
  131071,
  262139,
  524287,
  1048573,
  2097143,
  4194301,
  8388593,
  16777213,
  33554393,
  67108859,
  134217689,
  268435399,
  536870909,
  1073741789,
  2147483647  /* Maximum integer */
};

hecmw_hash_p *hecmw_hash_p_new(unsigned int index)
{
  unsigned int i, size;
  Bin *bin;
  hecmw_hash_p *hash;

  hash = (hecmw_hash_p *)malloc(sizeof(hecmw_hash_p));
  if (hash == NULL) return NULL;

  hash->num_put = 0;
  size = hecmw_hashsize_template[hecmw_hashsize_init_index];
  hash->n = size;
  hash->bin = (Bin *)malloc(hash->n*sizeof(Bin));
  if (hash->bin == NULL) {
    free(hash);
    return NULL;
  }

  bin = hash->bin;
  for (i=0; i< size; i++) {
    bin->n = 0;
    bin->list = NULL;
    bin++;
  }
  
  return hash;
}

void hecmw_hash_p_delete(hecmw_hash_p *hash)
{
  unsigned int i, j, k, l;
  Bin *bin;
  List *list;
  if (hash == NULL)  return;
  k = hash->n;
  bin = hash->bin;
  for (i=0; i< k; i++) {
    l = bin->n;
    list = bin->list;
    for (j=0; j<l; j++) {
      free(list->key);
      list++;
    }
    free(bin->list);
    bin++;
  }
  free(hash->bin);
  free(hash);
}

void *hecmw_hash_p_get(const hecmw_hash_p *hash, const char *key)
{
  unsigned int index;
  List *list;
  
  if (hash == NULL) return NULL;
  if (key == NULL) return NULL;

  index = hash_key(key) % hash->n;
  list = get_list(&(hash->bin[index]), key);

  if (list == NULL) return NULL;

  return list->value;
}

int hecmw_hash_p_exist(const hecmw_hash_p *hash, const char *key)
{
  unsigned int index;
  List *list;

  if (hash == NULL)  return 0;
  if (key == NULL)  return 0;

  index = hash_key(key) % hash->n;
  list = get_list(&(hash->bin[index]), key);

  if (list == NULL)  return 0;
  return 1;
}

int hecmw_hash_p_put(hecmw_hash_p *hash, const char *key, void *value)
{
  unsigned int key_len, index, stat;
  unsigned int hashkey, n;
  Bin *bin;
  List *tmp_list, *list;
  char *new_key;

  if (hash  == NULL) return 0;
  if (key   == NULL) return 0;
  if (value == NULL) return 0;
  
  stat = 0.8 * hash->n;
  if ( hash->num_put >= stat) {
    if (hecmw_hash_p_resize(hash) == 1) {
      return 1;
    }
  }
  hashkey = hash_key(key);
  index = hashkey % hash->n;
  key_len = strlen(key);

  bin = &(hash->bin[index]);
  list = get_list(bin, key);
  
  if (list != NULL) {
    return 1;
  }
  
  new_key = (char *)malloc((key_len+1)*sizeof(char));
  if (new_key == NULL) return 0;
  
  n = bin->n;
  if (n == 0) {
    tmp_list = (List *)malloc(sizeof(List));
    if (tmp_list == NULL) {
      free(new_key);
      return 0;
    }
    bin->list = tmp_list;
  } else {
    tmp_list = (List *)realloc(bin->list, (n+1)*sizeof(List));
    if (tmp_list == NULL) {
      free(new_key);
      return 0;
    }
    bin->list = tmp_list;
  }
  
  list = &(bin->list[n]);
  list->key = new_key;
  strcpy(list->key, key);
  list->value = value;
  list->hashkey = hashkey;
  bin->n++;
  hash->num_put++;
  
  return 1;
}

int hecmw_hash_p_resize(hecmw_hash_p *hash)
{
  unsigned int i, j, n;
  unsigned int newsize, index;
  Bin *newbin, *bin, *tmp_bin;
  List *newlist, *list, *tmp_list;
  
  j = 0;
  for (i=0; i<22; i++) {
    if (hash->n == hecmw_hashsize_template[i]) {
      j = i+1;
      break;
    }
  }

  if (j == 0) {
    printf("ERROR in hash table resize\n");
    return 1;
  }

  newsize = hecmw_hashsize_template[j];
  newbin = (Bin *)malloc(newsize*sizeof(Bin));
  if (newbin == NULL) return 1;

  tmp_bin = newbin;
  for (i=0; i< newsize; i++) {
    tmp_bin->n = 0;
    tmp_bin->list = NULL;
    tmp_bin++;
  }

  bin = hash->bin;
  for (i=0; i < hash->n; i++) { 
    list = bin->list;
    for (j=0; j<bin->n; j++) {
      index = list->hashkey % newsize;
      tmp_bin = &(newbin[index]);
      
      n = tmp_bin->n;
      if (n == 0) {
        newlist = (List *)malloc(sizeof(List));
        if (newlist == NULL) return 1;
        tmp_bin->list = newlist;
      } else {
        newlist = (List *)realloc(tmp_bin->list, (n+1)*sizeof(List));
        if (newlist == NULL) return 1;
        tmp_bin->list = newlist;
      }
      
      tmp_list = &(tmp_bin->list[n]);
      tmp_list->key = list->key;
      tmp_list->value = list->value;
      tmp_list->hashkey = list->hashkey;
      tmp_bin->n++;
      list++;
    }
    free(bin->list);
    bin++;
  }
  
  free(hash->bin);
  hash->bin = newbin;
  hash->n = newsize;

  return 0;
}

static List *get_list(Bin *bin, const char *key)
{
  unsigned int i, n;
  List *list;

  n = bin->n;
  if (n == 0) return NULL;

  list = bin->list;
  for (i=0; i<n; i++) {
    if (list->key != NULL && list->value != NULL) {
      if (strcmp(list->key, key) == 0) {
        return list;
      }
    }
    list++;
  }
  return NULL;
}

/** djb2 hash function **/
static unsigned int hash_key(const char *c)
{
  unsigned int hash_key = 5381;
  int i;

  while (i = *c++) {
    hash_key = ((hash_key << 5) + hash_key) + i;
  }

  return hash_key;
}
