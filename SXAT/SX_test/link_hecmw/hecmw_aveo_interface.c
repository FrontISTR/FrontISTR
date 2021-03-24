#include <ve_offload.h>
#include <stdlib.h>
#include <stdio.h>

void check_ve(uint64_t ret){
    if(ret != 0){
        printf("Error! ret = %d\n", ret);
//        abort();
    }
}

extern int get_aveo_version()
{
    struct veo_proc_handle *proc = veo_proc_create_static(0, "../../hecmw_ve/hecmw_ve");
    uint64_t handle = NULL;/* find a function in the executable */
    struct veo_thr_ctxt *ctx = veo_context_open(proc);
    struct veo_args *argp = veo_args_alloc();
    uint64_t retval;

    int ret;
    size_t N = 100;
    size_t val = 100;

    // vector x
    double* x = (double*)malloc(sizeof(double)*N);
    for(int i=0; i<N; i++){
        x[i] = i;
    }

    uint64_t xd;
    ret = veo_alloc_mem(proc, &xd, N * sizeof(double)); check_ve(ret);
    ret = veo_write_mem(proc, xd, x, N * sizeof(double)); check_ve(ret);
//   ret = veo_read_mem(proc, x, xd, N * sizeof(double));

    long y = 0;
    ret = veo_args_set_u64(argp, 0, val); check_ve(ret);
    ret = veo_args_set_u64(argp, 1, N); check_ve(ret);
    ret = veo_args_set_u64(argp, 2, xd); check_ve(ret);
    printf("goma(p): %p\n", xd);
    printf("goma(d): %ld\n", xd);
//     veo_args_set_stack(argp, VEO_INTENT_IN, 0, (char*)&xd, sizeof(xd));
//     veo_args_set_stack(argp, VEO_INTENT_IN, 1, (char*)&N, sizeof(N));
//     veo_args_set_stack(argp, VEO_INTENT_OUT, 2, (char*)&y, sizeof(y));
// 
// 
    const char func[] = "hecmw_solver_iterative_sxat_MP_hecmw_solve_sxat_entry_";
    uint64_t id = veo_call_async_by_name(ctx, handle, func, argp);

    uint64_t req;
    ret = veo_call_wait_result(ctx, req, &retval);
    //check_ve(ret);

//   ret = veo_read_mem(proc, y, yd, N * sizeof(double));

    //////////

//     printf("===VH====\n");
//     for(size_t i = 0; i < N; i++){
//         printf("%f\n",y[i]);
//     }

    veo_args_free(argp);
    veo_context_close(ctx);
    veo_proc_destroy(proc);
    return 123;
}
