/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#undef NDEBUG
#include <stdio.h>
#include <assert.h>
#include "api/hecmw_api_result_data.h"

int main()
{
    printf("FrontISTR API Test : hecmwST_result_data\n");
    void* result = hecmw_api_result_new();
    assert(result != 0);
    hecmw_api_result_delete(result);
    return 0;
}