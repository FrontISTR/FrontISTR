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

#include <stdlib.h>
#include <string.h>
#include "hecmw_io_struct.h"

typedef struct hecmw_hash_p hecmw_hash_p;

hecmw_hash_p * hecmw_hash_p_new(unsigned int capacity);

void hecmw_hash_p_delete(hecmw_hash_p *hash);

void *hecmw_hash_p_get(const hecmw_hash_p *hash, const char *key);

int hecmw_hash_p_exist(const hecmw_hash_p *hash, const char *key);

int hecmw_hash_p_put(hecmw_hash_p *hash, const char *key, void *value);


