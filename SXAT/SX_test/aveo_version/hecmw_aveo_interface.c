#include <ve_offload.h>
#include <stdlib.h>
#include <stdio.h>

extern int get_aveo_version()
{
    printf("%d\n",veo_api_version());
    return veo_api_version();
}
