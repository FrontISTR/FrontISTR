#include<stdio.h>
#include<stdlib.h>

#if defined USE_SX
#include<ve_offload.h>
#endif

extern void hecmw_solve_sxat_test_fortran();

int hecmw_solve_sxat_test(double value, size_t size, double* array){

    printf("[VE] val = %f, size = %d\n", value, size);

    size_t i;
    for(i=0; i<size; i++){
        array[i] += value;
    }

    hecmw_solve_sxat_test_fortran(value, size, array);
    return 0;
}

#if defined USE_SX
void check_ve(uint64_t ret){
    if(ret != 0){
        printf("Error! ret = %d\n", ret);
//        abort();
    }
}
#endif

void hecmw_solve_cg_entry(
        int N, //0
        int NP, //1
        int NDOF, //2
        double* D, //3
        double* AL, //4
        double* AU, //5
        int* indexL, //6
        int* indexU, //7
        int* itemL, //8
        int* itemU, //9
        int my_rank, //10
        double* x, //11
        double* b, //12
        int iterlog, //13
        int timelog, //14
        int estcond, //15
        int iter, //16
        int error, //17
        double Tset, //18
        double Tsol, //19
        double Tcomm //20
        ){
#if defined USE_SX
    //printf("[VH] I will call ve!\n");
    const char* veobj = getenv("FRONTISTR_VEOBJ");
    struct veo_proc_handle *proc = veo_proc_create_static(0, veobj);
    uint64_t handle = NULL;/* find a function in the executable */
    struct veo_thr_ctxt *ctx = veo_context_open(proc);
    struct veo_args *argp = veo_args_alloc();

    int ret;
    double val = 123;
    int i;

    int size = N * NDOF;
    int nnz = N * NDOF * NP;

    uint64_t Dd;
    veo_alloc_mem(proc, &Dd, size * sizeof(double)); check_ve(ret);
    veo_write_mem(proc, Dd, D, size * sizeof(double)); check_ve(ret);
    uint64_t ALd;
    veo_alloc_mem(proc, &ALd, nnz * sizeof(double)); check_ve(ret);
    veo_write_mem(proc, ALd, AL, nnz * sizeof(double)); check_ve(ret);
    uint64_t AUd;
    veo_alloc_mem(proc, &AUd, nnz * sizeof(double)); check_ve(ret);
    veo_write_mem(proc, AUd, AU, nnz * sizeof(double)); check_ve(ret);
    uint64_t indexLd;
    veo_alloc_mem(proc, &indexLd, nnz * sizeof(int)); check_ve(ret);
    veo_write_mem(proc, indexLd, indexL, size * sizeof(int)); check_ve(ret);
    uint64_t indexUd;
    veo_alloc_mem(proc, &indexUd, nnz * sizeof(int)); check_ve(ret);
    veo_write_mem(proc, indexUd, indexU, size * sizeof(int)); check_ve(ret);
    uint64_t itemLd;
    veo_alloc_mem(proc, &itemLd, size * sizeof(int)); check_ve(ret);
    veo_write_mem(proc, itemLd, itemL, size * sizeof(int)); check_ve(ret);
    uint64_t itemUd;
    veo_alloc_mem(proc, &itemUd, size * sizeof(int)); check_ve(ret);
    veo_write_mem(proc, itemUd, itemU, size * sizeof(int)); check_ve(ret);
    uint64_t xd;
    veo_alloc_mem(proc, &xd, size * sizeof(double)); check_ve(ret);
    veo_write_mem(proc, xd, x, size * sizeof(double)); check_ve(ret);
    uint64_t xd;
    veo_alloc_mem(proc, &xd, size * sizeof(double)); check_ve(ret);
    veo_write_mem(proc, xd, x, size * sizeof(double)); check_ve(ret);
    uint64_t bd;
    veo_alloc_mem(proc, &bd, size * sizeof(double)); check_ve(ret);
    veo_write_mem(proc, bd, b, size * sizeof(double)); check_ve(ret);


    ret = veo_args_set_i32(argp, 0, N); check_ve(ret);
    ret = veo_args_set_i32(argp, 1, NP); check_ve(ret);
    ret = veo_args_set_i32(argp, 2, NDOF); check_ve(ret);
    ret = veo_args_set_u64(argp, 3, Dd); check_ve(ret);
    ret = veo_args_set_u64(argp, 4, ALd); check_ve(ret);
    ret = veo_args_set_u64(argp, 5, AUd); check_ve(ret);
    ret = veo_args_set_u64(argp, 6, indexLd); check_ve(ret);
    ret = veo_args_set_u64(argp, 7, indexUd); check_ve(ret);
    ret = veo_args_set_u64(argp, 8, itemLd); check_ve(ret);
    ret = veo_args_set_u64(argp, 9, itemUd); check_ve(ret);
    ret = veo_args_set_i32(argp, 10, my_rank); check_ve(ret);
    ret = veo_args_set_u64(argp, 11, xd); check_ve(ret);
    ret = veo_args_set_u64(argp, 12, bd); check_ve(ret);
    ret = veo_args_set_i32(argp, 13, iterlog); check_ve(ret);
    ret = veo_args_set_i32(argp, 14, timelog); check_ve(ret);
    ret = veo_args_set_i32(argp, 15, estcond); check_ve(ret);
    ret = veo_args_set_i32(argp, 16, iter); check_ve(ret);
    ret = veo_args_set_i32(argp, 17, error); check_ve(ret);
    ret = veo_args_set_double(argp, 18, Tset); check_ve(ret);
    ret = veo_args_set_double(argp, 19, Tsol); check_ve(ret);
    ret = veo_args_set_double(argp, 20, Tcomm); check_ve(ret);

    const char func[] = "hecmw_solve_sxat_cg_interface";
    uint64_t id = veo_call_async_by_name(ctx, handle, func, argp);

    uint64_t req;
    uint64_t retval;
    ret = veo_call_wait_result(ctx, req, &retval);
    //check_ve(ret);
    //printf("ret is %d\n",retval);

   ret = veo_read_mem(proc, x, xd, N * sizeof(double));

//     printf("==========VH==========\n");
//     for(i = 0; i < N; i++){
//         printf("%f\n",x[i]);
//     }

    veo_args_free(argp);
    veo_context_close(ctx);
    veo_proc_destroy(proc);
#endif
}

