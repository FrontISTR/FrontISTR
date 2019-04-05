/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "hecmw_io_struct.h"

typedef struct hecmw_hash_p hecmw_hash_p;

hecmw_hash_p *hecmw_hash_p_new(unsigned int capacity);

void hecmw_hash_p_delete(hecmw_hash_p *hash);

void *hecmw_hash_p_get(const hecmw_hash_p *hash, const char *key);

int hecmw_hash_p_exist(const hecmw_hash_p *hash, const char *key);

int hecmw_hash_p_put(hecmw_hash_p *hash, const char *key, void *value);
