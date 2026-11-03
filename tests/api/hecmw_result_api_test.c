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
    void* result1 = hecmw_api_result_new();
    void* result2 = hecmw_api_result_new();
    assert(result1 != 0);
    assert(result2 != 0);
    assert(hecmw_api_result_ng_component(result1)==0);
    assert(hecmw_api_result_nn_component(result1)==0);
    assert(hecmw_api_result_ne_component(result1)==0);
    hecmw_api_result_delete(result1);
    assert(hecmw_api_result_ng_component(result2)==0);
    assert(hecmw_api_result_nn_component(result2)==0);
    assert(hecmw_api_result_ne_component(result2)==0);
    hecmw_api_result_delete(result2);
    return 0;
}
