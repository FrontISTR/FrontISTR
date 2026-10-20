/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

 #include <stdio.h>
#include "api/fstr_api_param.h"

int main(void)
{
    printf("FrontISTR API Test : fstr_param\n");
    void* param = fstr_api_param_new();
    fstr_api_param_delete(param);
    return 0;
}