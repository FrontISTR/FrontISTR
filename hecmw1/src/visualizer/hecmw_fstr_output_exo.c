/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *
 * Exodus II output via NetCDF API (without libexodus dependency)
 *****************************************************************************/
#include "hecmw_fstr_output_exo.h"

#ifdef WITH_NETCDF

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <netcdf.h>
#include "hecmw_malloc.h"
#include "hecmw_etype.h"
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_comm_util.h"

/*============================================================================
 * Constants
 *============================================================================*/
#define EXO_MAX_STR_LENGTH  33
#define EXO_MAX_LINE_LENGTH 81
#define EXO_API_VERSION     8.03f
#define EXO_DB_VERSION      8.03f
#define EXO_FILE_SIZE       1

/* Maximum number of distinct element types we handle */
#define MAX_ELEM_BLOCKS 64

/*============================================================================
 * FrontISTR element type -> Exodus element type string mapping
 *
 * FrontISTR uses numeric element type codes defined in hecmw_common_define.h.
 * Exodus requires a string like "TETRA4".
 *
 * The actual number of nodes written to connectivity is determined at runtime
 * from elem_node_index (with shift applied), matching the VTK output behavior
 * exactly. The exo_name here tells ParaView how to interpret the topology.
 *
 * Element categories:
 *   1xx: Rod/Bar       2xx: 2D surface    3xx: 3D solid
 *   4xx: Master/slave  5xx: Joint         6xx: Beam
 *   7xx: Shell         9xx: Interface(LN) 10xx: Patch(PT)
 *============================================================================*/
typedef struct {
    int    hecmw_type;      /* FrontISTR element type code */
    const char *exo_name;   /* Exodus element type string  */
} ElemTypeMap;

static const ElemTypeMap elem_type_table[] = {
    /* --- 1D Rod/Bar elements --- */
    {111,   "BAR2"},       /* ROD1:  2 nodes */
    {112,   "BAR3"},       /* ROD2:  3 nodes */
    {301,   "BAR2"},       /* ROD31: 2 nodes */
    /* --- 2D Surface elements --- */
    {231,   "TRI3"},       /* TRI1:  3 nodes */
    {232,   "TRI6"},       /* TRI2:  6 nodes */
    {2322,  "TRI3"},       /* TRI22: variant, VTK treats as TRI */
    {241,   "QUAD4"},      /* QUA1:  4 nodes */
    {242,   "QUAD8"},      /* QUA2:  8 nodes */
    /* --- 3D Solid elements --- */
    {341,   "TETRA4"},     /* TET1:   4 nodes */
    {3414,  "TETRA4"},     /* TET1_4: 4 nodes (variant) */
    {342,   "TETRA10"},    /* TET2:  10 nodes (node reorder needed) */
    {3422,  "TETRA4"},     /* TET22:  VTK treats as TETRA (4 nodes) */
    {351,   "WEDGE6"},     /* PRI1:   6 nodes */
    {352,   "WEDGE15"},    /* PRI2:  15 nodes */
    {361,   "HEX8"},       /* HEX1:   8 nodes */
    {3614,  "HEX8"},       /* HEX1_4: 8 nodes (variant) */
    {362,   "HEX20"},      /* HEX2:  20 nodes */
    {371,   "PYRAMID5"},   /* PYR1:   5 nodes */
    {372,   "PYRAMID13"},  /* PYR2:  13 nodes */
    /* --- Master/Slave surface elements --- */
    {431,   "TETRA4"},     /* MST1: 4 nodes, VTK=TRI but stored as tet */
    {432,   "TETRA4"},     /* MST2: 7 nodes, VTK=TRI */
    {441,   "PYRAMID5"},   /* MSQ1: 5 nodes, VTK=QUAD but stored as pyr */
    {442,   "PYRAMID5"},   /* MSQ2: 9 nodes, VTK=QUAD */
    /* --- Joint elements --- */
    {501,   "BAR2"},       /* JTB1:    2 nodes */
    {511,   "BAR2"},       /* SPGDPT1: 2 nodes */
    {531,   "TRI3"},       /* JTT1:  6 nodes total, VTK=TRI (uses first 3) */
    {532,   "TRI3"},       /* JTT2: 12 nodes total, VTK=TRI */
    {541,   "QUAD4"},      /* JTQ1:  8 nodes total, VTK=QUAD (uses first 4) */
    {542,   "QUAD4"},      /* JTQ2: 16 nodes total, VTK=QUAD */
    /* --- Beam elements --- */
    {611,   "BAR2"},       /* BEM1: 2 nodes */
    {612,   "BAR3"},       /* BEM2: 3 nodes */
    {641,   "BAR2"},       /* BEM3: 4 nodes, shift=2 -> 2 nodes (mixed beam-341) */
    /* --- Shell elements --- */
    {731,   "TRI3"},       /* SHT1: 3 nodes */
    {732,   "TRI6"},       /* SHT2: 6 nodes */
    {741,   "QUAD4"},      /* SHQ1: 4 nodes */
    {742,   "QUAD8"},      /* SHQ2: 8 nodes */
    {743,   "QUAD9"},      /* SHQ3: 9 nodes */
    {761,   "TRI3"},       /* SHT6: 6 nodes, shift=3 -> 3 nodes (mixed shell-solid) */
    {781,   "QUAD4"},      /* SHQ8: 8 nodes, shift=4 -> 4 nodes (mixed shell-solid) */
    /* --- Interface/Link elements (LN) --- all 2-node lines */
    {911,   "BAR2"}, {912, "BAR2"}, {913, "BAR2"}, {914, "BAR2"}, {915, "BAR2"}, {916, "BAR2"},
    {921,   "BAR2"}, {922, "BAR2"}, {923, "BAR2"}, {924, "BAR2"}, {925, "BAR2"}, {926, "BAR2"},
    {931,   "BAR2"}, {932, "BAR2"}, {933, "BAR2"}, {934, "BAR2"}, {935, "BAR2"}, {936, "BAR2"},
    {941,   "BAR2"}, {942, "BAR2"}, {943, "BAR2"}, {944, "BAR2"}, {945, "BAR2"}, {946, "BAR2"},
    {951,   "BAR2"}, {952, "BAR2"}, {953, "BAR2"}, {954, "BAR2"}, {955, "BAR2"}, {956, "BAR2"},
    {961,   "BAR2"}, {962, "BAR2"}, {963, "BAR2"}, {964, "BAR2"}, {965, "BAR2"}, {966, "BAR2"},
    /* --- Patch elements (PT) --- */
    {1031,  "TRI3"},       /* PTT1: 3 nodes */
    {1032,  "TRI6"},       /* PTT2: 6 nodes */
    {1041,  "QUAD4"},      /* PTQ1: 4 nodes */
    {1042,  "QUAD8"},      /* PTQ2: 8 nodes */
    /* --- Sentinel --- */
    {0,     NULL}
};

/*============================================================================
 * Helper: NetCDF error check macro
 *============================================================================*/
#define NCERR(e) do { \
    if ((e) != NC_NOERR) { \
        fprintf(stderr, "NetCDF error at %s:%d: %s\n", \
                __FILE__, __LINE__, nc_strerror(e)); \
        return; \
    } \
} while(0)

/*============================================================================
 * Helper: look up Exodus element type string for a FrontISTR type code
 *============================================================================*/
static const ElemTypeMap* get_elem_type_map(int hecmw_type)
{
    int i;
    for (i = 0; elem_type_table[i].exo_name != NULL; i++) {
        if (elem_type_table[i].hecmw_type == hecmw_type)
            return &elem_type_table[i];
    }
    return NULL;
}

/*============================================================================
 * Helper: get node shift for special element types (same as VTK output)
 *============================================================================*/
static int get_node_shift(int hecmw_type)
{
    switch (hecmw_type) {
        case 641: return 2;
        case 761: return 3;
        case 781: return 4;
        default:  return 0;
    }
}

/*============================================================================
 * Helper: get the effective number of nodes for an element
 *         (after applying shift, matching VTK output behavior)
 *============================================================================*/
static int get_effective_n_nodes(struct hecmwST_local_mesh *mesh, int elem_idx)
{
    int jS = mesh->elem_node_index[elem_idx];
    int jE = mesh->elem_node_index[elem_idx + 1];
    int shift = get_node_shift(mesh->elem_type[elem_idx]);
    return jE - jS - shift;
}

/*============================================================================
 * Element block info (gathered by scanning mesh element types)
 *============================================================================*/
typedef struct {
    int  hecmw_type;           /* FrontISTR element type code */
    const char *exo_name;      /* Exodus element type string */
    int  nod_per_elem;         /* nodes per element */
    int  num_elem;             /* number of elements in this block */
    int *elem_indices;         /* indices into original element array */
} ElemBlockInfo;

/*============================================================================
 * Helper: write a padded string into a char array for NetCDF
 *============================================================================*/
static void pad_string(char *buf, int len, const char *str)
{
    int slen = (int)strlen(str);
    int i;
    if (slen > len) slen = len;
    for (i = 0; i < slen; i++) buf[i] = str[i];
    for (i = slen; i < len; i++) buf[i] = '\0';
}

/*============================================================================
 * Helper: expand multi-component FrontISTR variable labels into per-component
 *         Exodus variable names.
 *
 * FrontISTR stores a single label like "DISPLACEMENT" with nn_dof[i]=3,
 * meaning 3 components (X,Y,Z). Exodus needs individual names per component:
 *   "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"
 *
 * For 6-component tensor data (e.g., stress):
 *   "STRESS_XX", "STRESS_YY", "STRESS_ZZ", "STRESS_XY", "STRESS_YZ", "STRESS_ZX"
 *============================================================================*/
static const char *comp_suffix_1[] = {""};
static const char *comp_suffix_2[] = {"_X", "_Y"};
static const char *comp_suffix_3[] = {"_X", "_Y", "_Z"};
static const char *comp_suffix_3_principal[] = {"_1st", "_2nd", "_3rd"};
static const char *comp_suffix_6[] = {"_XX", "_YY", "_ZZ", "_XY", "_YZ", "_ZX"};
static const char *comp_suffix_7[] = {"_XX", "_YY", "_ZZ", "_XY", "_YZ", "_ZX", "_MISES"};

static int label_is_principal(const char *label)
{
    /* Match FrontISTR labels:
     *   NodalPrincipalSTRESS, NodalPrincipalSTRAIN,
     *   ElementalPrincipalSTRESS, ElementalPrincipalSTRAIN, etc. */
    return (strstr(label, "Principal") != NULL);
}

static const char** get_comp_suffixes(int ndof, const char *label)
{
    switch (ndof) {
        case 1: return comp_suffix_1;
        case 2: return comp_suffix_2;
        case 3: return label_is_principal(label) ? comp_suffix_3_principal
                                                 : comp_suffix_3;
        case 6: return comp_suffix_6;
        case 7: return comp_suffix_7;
        default: return NULL;
    }
}

/*============================================================================
 * Main Exodus output function
 *
 * This function writes an Exodus II file using NetCDF API calls directly.
 * It follows the same calling convention as vtk_output / bin_vtk_output.
 *
 * On the first call (time_step_count == 0), it creates the file and defines
 * the full schema. On subsequent calls, it opens the file and appends a
 * new time step.
 *
 * File naming:
 *   - Each rank writes: {outfile1}/{outfile}.{myrank}.exo
 *   (No pvtu-equivalent needed: ParaView can load .exo.N.M Nemesis files,
 *    or individual per-rank .exo files can be loaded as a group.)
 *============================================================================*/
void exodus_output(struct hecmwST_local_mesh *mesh,
                   struct hecmwST_result_data *data,
                   char *outfile, char *outfile1,
                   HECMW_Comm VIS_COMM,
                   int per_step)
{
    int i, j, k;
    int jS, jE;
    int myrank, petot;
    int n_node, n_elem, shift;
    int data_tot_n, data_tot_e;
    int table342[10] = {0, 1, 2, 3, 6, 4, 5, 7, 8, 9};
    char *p;

    int ncid, stat;
    int old_fill_mode;
    static int time_step_count = 0;

    /*
     * Two modes:
     *   per_step == 0 (EXODUS):      All timesteps in one file (append mode)
     *   per_step == 1 (STEP_EXODUS): One file per timestep (always create new)
     *
     * For per_step==0, the file path is fixed across calls (static).
     * For per_step==1, a new file is created each call with the step number.
     */
    static char file_exo[HECMW_FILENAME_LEN];
    static char dir_exo[HECMW_FILENAME_LEN];
    static char basename_exo[HECMW_FILENAME_LEN];

    /* --- Element block working data --- */
    int num_blocks = 0;
    ElemBlockInfo blocks[MAX_ELEM_BLOCKS];

    /* --- Count total DOFs for nodal/element data --- */
    int total_nod_var = 0;  /* total number of scalar Exodus nodal variables */
    int total_elem_var = 0; /* total number of scalar Exodus element variables */

    HECMW_Comm_rank(VIS_COMM, &myrank);
    HECMW_Comm_size(VIS_COMM, &petot);
    n_node = mesh->n_node;
    n_elem = mesh->n_elem;

    data_tot_n = 0;
    for (i = 0; i < data->nn_component; i++)
        data_tot_n += data->nn_dof[i];

    data_tot_e = 0;
    for (i = 0; i < data->ne_component; i++)
        data_tot_e += data->ne_dof[i];

    /* Count total scalar variables for Exodus (expand multi-component) */
    for (i = 0; i < data->nn_component; i++)
        total_nod_var += data->nn_dof[i];
    for (i = 0; i < data->ne_component; i++)
        total_elem_var += data->ne_dof[i];

    /* --- Build file path --- */
    if (time_step_count == 0) {
        /*
         * First call: compute dir_exo and basename_exo from caller's args.
         *
         * outfile1 is e.g. "vis_out_psf.0001" (with timestep suffix).
         * We strip the suffix to get a stable directory: "vis_out_psf_exo".
         */

        /* Build directory: strip timestep suffix from outfile1 and append _exo */
        strncpy(dir_exo, outfile1, HECMW_FILENAME_LEN - 1);
        dir_exo[HECMW_FILENAME_LEN - 1] = '\0';
        p = strrchr(dir_exo, '.');
        if (p != NULL) *p = '\0'; /* "vis_out_psf.0001" -> "vis_out_psf" */
        strncat(dir_exo, "_exo", HECMW_FILENAME_LEN - strlen(dir_exo) - 1);
        /* dir_exo = "vis_out_psf_exo" */

        /* Build basename: strip directory prefix and timestep suffix from outfile */
        strncpy(basename_exo, outfile, HECMW_FILENAME_LEN - 1);
        basename_exo[HECMW_FILENAME_LEN - 1] = '\0';
        p = strrchr(basename_exo, '/');
        if (p != NULL) {
            memmove(basename_exo, p + 1, strlen(p + 1) + 1);
        }
        p = strrchr(basename_exo, '.');
        if (p != NULL) *p = '\0'; /* "vis_out_psf.0001" -> "vis_out_psf" */
    }

    if (per_step) {
        /*
         * STEP_EXODUS mode: one file per timestep
         *   serial:   {basename}.step{NNNN}.exo
         *   parallel: {basename}.step{NNNN}.exo.{nprocs}.{rank}  (Nemesis)
         */
        if (petot == 1) {
            snprintf(file_exo, HECMW_FILENAME_LEN, "%s/%s.step%04d.exo",
                    dir_exo, basename_exo, time_step_count);
        } else {
            snprintf(file_exo, HECMW_FILENAME_LEN, "%s/%s.step%04d.exo.%d.%d",
                    dir_exo, basename_exo, time_step_count, petot, myrank);
        }
    } else {
        /*
         * EXODUS mode: all timesteps in one file (only set path on first call)
         *   serial:   {basename}.exo
         *   parallel: {basename}.exo.{nprocs}.{rank}  (Nemesis)
         */
        if (time_step_count == 0) {
            if (petot == 1) {
                snprintf(file_exo, HECMW_FILENAME_LEN, "%s/%s.exo",
                        dir_exo, basename_exo);
            } else {
                snprintf(file_exo, HECMW_FILENAME_LEN, "%s/%s.exo.%d.%d",
                        dir_exo, basename_exo, petot, myrank);
            }
        }
        /* For subsequent calls, file_exo is already set (static) */
    }

    /* Create the directory (needed on first call, or every call for per_step) */
    if (time_step_count == 0 || per_step) {
        if (HECMW_ctrl_make_subdir(file_exo)) {
            HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0009: Cannot open output directory for Exodus");
        }
    }

    /*========================================================================
     * PASS 1: Scan elements and group into blocks by element type
     *========================================================================*/
    memset(blocks, 0, sizeof(blocks));
    for (i = 0; i < n_elem; i++) {
        int etype = mesh->elem_type[i];
        int found = 0;
        for (j = 0; j < num_blocks; j++) {
            if (blocks[j].hecmw_type == etype) {
                blocks[j].num_elem++;
                found = 1;
                break;
            }
        }
        if (!found) {
            if (num_blocks >= MAX_ELEM_BLOCKS) {
                HECMW_vis_print_exit("ERROR: Too many element types for Exodus output");
            }
            const ElemTypeMap *map = get_elem_type_map(etype);
            if (map == NULL) {
                fprintf(stderr, "WARNING: Unknown element type %d, skipping\n", etype);
                continue;
            }
            blocks[num_blocks].hecmw_type = etype;
            blocks[num_blocks].exo_name = map->exo_name;
            /* Determine nodes per element from first element of this type */
            blocks[num_blocks].nod_per_elem = get_effective_n_nodes(mesh, i);
            blocks[num_blocks].num_elem = 1;
            num_blocks++;
        }
    }

    /* Allocate index arrays and fill element indices for each block */
    for (j = 0; j < num_blocks; j++) {
        blocks[j].elem_indices = (int *)HECMW_malloc(sizeof(int) * blocks[j].num_elem);
        blocks[j].num_elem = 0; /* reset count, will re-fill */
    }
    for (i = 0; i < n_elem; i++) {
        int etype = mesh->elem_type[i];
        for (j = 0; j < num_blocks; j++) {
            if (blocks[j].hecmw_type == etype) {
                blocks[j].elem_indices[blocks[j].num_elem] = i;
                blocks[j].num_elem++;
                break;
            }
        }
    }

    /*========================================================================
     * CREATE or OPEN file
     *
     * per_step==1: Always create a new file (one file per timestep)
     * per_step==0: Create on first call, open for append on subsequent calls
     *========================================================================*/
    if (per_step || time_step_count == 0) {
        /*====================================================================
         * FIRST CALL: Create file and define full schema
         *====================================================================*/
        stat = nc_create(file_exo, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
        NCERR(stat);
        stat = nc_set_fill(ncid, NC_NOFILL, &old_fill_mode);
        NCERR(stat);

        /*----------------------------------------------------------------
         * Global Attributes  [REQUIRED]
         *----------------------------------------------------------------*/
        {
            float api_ver = EXO_API_VERSION;
            float db_ver  = EXO_DB_VERSION;
            int   word_sz = (int)sizeof(double);
            int   file_sz = EXO_FILE_SIZE;
            char  title[] = "FrontISTR Exodus output";

            stat = nc_put_att_float(ncid, NC_GLOBAL, "api_version", NC_FLOAT, 1, &api_ver);
            NCERR(stat);
            stat = nc_put_att_float(ncid, NC_GLOBAL, "version", NC_FLOAT, 1, &db_ver);
            NCERR(stat);
            stat = nc_put_att_int(ncid, NC_GLOBAL, "floating_point_word_size", NC_INT, 1, &word_sz);
            NCERR(stat);
            stat = nc_put_att_int(ncid, NC_GLOBAL, "file_size", NC_INT, 1, &file_sz);
            NCERR(stat);
            stat = nc_put_att_text(ncid, NC_GLOBAL, "title", strlen(title), title);
            NCERR(stat);
        }

        /*----------------------------------------------------------------
         * Dimensions  [REQUIRED]
         *----------------------------------------------------------------*/
        int dim_len_string, dim_len_name, dim_len_line, dim_four;
        int dim_time_step, dim_num_dim, dim_num_nodes, dim_num_elem;
        int dim_num_el_blk;
        int dim_num_qa_rec;

        stat = nc_def_dim(ncid, "len_string", EXO_MAX_STR_LENGTH, &dim_len_string);
        NCERR(stat);
        stat = nc_def_dim(ncid, "len_name",   EXO_MAX_STR_LENGTH, &dim_len_name);
        NCERR(stat);
        stat = nc_def_dim(ncid, "len_line",   EXO_MAX_LINE_LENGTH, &dim_len_line);
        NCERR(stat);
        stat = nc_def_dim(ncid, "four",       4, &dim_four);
        NCERR(stat);
        stat = nc_def_dim(ncid, "time_step",  NC_UNLIMITED, &dim_time_step);
        NCERR(stat);
        stat = nc_def_dim(ncid, "num_dim",    3, &dim_num_dim);
        NCERR(stat);
        stat = nc_def_dim(ncid, "num_nodes",  (size_t)n_node, &dim_num_nodes);
        NCERR(stat);
        stat = nc_def_dim(ncid, "num_elem",   (size_t)n_elem, &dim_num_elem);
        NCERR(stat);
        stat = nc_def_dim(ncid, "num_el_blk", (size_t)num_blocks, &dim_num_el_blk);
        NCERR(stat);

        /* [OPTIONAL] QA record dimension */
        stat = nc_def_dim(ncid, "num_qa_rec", 1, &dim_num_qa_rec);
        NCERR(stat);

        /* Per-block dimensions */
        int dim_el_in_blk[MAX_ELEM_BLOCKS];
        int dim_nod_per_el[MAX_ELEM_BLOCKS];
        for (j = 0; j < num_blocks; j++) {
            char dname[EXO_MAX_STR_LENGTH];
            snprintf(dname, EXO_MAX_STR_LENGTH, "num_el_in_blk%d", j + 1);
            stat = nc_def_dim(ncid, dname, (size_t)blocks[j].num_elem, &dim_el_in_blk[j]);
            NCERR(stat);
            snprintf(dname, EXO_MAX_STR_LENGTH, "num_nod_per_el%d", j + 1);
            stat = nc_def_dim(ncid, dname, (size_t)blocks[j].nod_per_elem, &dim_nod_per_el[j]);
            NCERR(stat);
        }

        /* Variable count dimensions */
        int dim_num_nod_var = -1, dim_num_elem_var = -1;
        if (total_nod_var > 0) {
            stat = nc_def_dim(ncid, "num_nod_var", (size_t)total_nod_var, &dim_num_nod_var);
            NCERR(stat);
        }
        if (total_elem_var > 0) {
            stat = nc_def_dim(ncid, "num_elem_var", (size_t)total_elem_var, &dim_num_elem_var);
            NCERR(stat);
        }

        /*----------------------------------------------------------------
         * Variables: Time  [REQUIRED]
         *----------------------------------------------------------------*/
        int var_time_whole;
        stat = nc_def_var(ncid, "time_whole", NC_DOUBLE, 1, &dim_time_step, &var_time_whole);
        NCERR(stat);

        /*----------------------------------------------------------------
         * Variables: Coordinates  [REQUIRED]
         *----------------------------------------------------------------*/
        int var_coordx, var_coordy, var_coordz;
        stat = nc_def_var(ncid, "coordx", NC_DOUBLE, 1, &dim_num_nodes, &var_coordx);
        NCERR(stat);
        stat = nc_def_var(ncid, "coordy", NC_DOUBLE, 1, &dim_num_nodes, &var_coordy);
        NCERR(stat);
        stat = nc_def_var(ncid, "coordz", NC_DOUBLE, 1, &dim_num_nodes, &var_coordz);
        NCERR(stat);

        /* coor_names [REQUIRED] */
        int var_coor_names;
        {
            int dims2[2] = {dim_num_dim, dim_len_name};
            stat = nc_def_var(ncid, "coor_names", NC_CHAR, 2, dims2, &var_coor_names);
            NCERR(stat);
        }

        /*----------------------------------------------------------------
         * Variables: Element Block info  [REQUIRED]
         *----------------------------------------------------------------*/
        int var_eb_status, var_eb_prop1;
        stat = nc_def_var(ncid, "eb_status", NC_INT, 1, &dim_num_el_blk, &var_eb_status);
        NCERR(stat);
        stat = nc_def_var(ncid, "eb_prop1",  NC_INT, 1, &dim_num_el_blk, &var_eb_prop1);
        NCERR(stat);
        stat = nc_put_att_text(ncid, var_eb_prop1, "name", 2, "ID");
        NCERR(stat);

        /* [OPTIONAL] eb_names */
        int var_eb_names;
        {
            int dims2[2] = {dim_num_el_blk, dim_len_name};
            stat = nc_def_var(ncid, "eb_names", NC_CHAR, 2, dims2, &var_eb_names);
            NCERR(stat);
        }

        /* Connectivity per block [REQUIRED] */
        int var_connect[MAX_ELEM_BLOCKS];
        for (j = 0; j < num_blocks; j++) {
            char vname[EXO_MAX_STR_LENGTH];
            snprintf(vname, EXO_MAX_STR_LENGTH, "connect%d", j + 1);
            int dims2[2] = {dim_el_in_blk[j], dim_nod_per_el[j]};
            stat = nc_def_var(ncid, vname, NC_INT, 2, dims2, &var_connect[j]);
            NCERR(stat);
            stat = nc_put_att_text(ncid, var_connect[j], "elem_type",
                                   strlen(blocks[j].exo_name), blocks[j].exo_name);
            NCERR(stat);
        }

        /*----------------------------------------------------------------
         * Variables: QA records  [OPTIONAL - nice to have]
         *----------------------------------------------------------------*/
        int var_qa_records;
        {
            int dims3[3] = {dim_num_qa_rec, dim_four, dim_len_string};
            stat = nc_def_var(ncid, "qa_records", NC_CHAR, 3, dims3, &var_qa_records);
            NCERR(stat);
        }

        /*----------------------------------------------------------------
         * Variables: ID maps  [OPTIONAL - nice to have for parallel]
         *----------------------------------------------------------------*/
        int var_node_num_map, var_elem_num_map;
        stat = nc_def_var(ncid, "node_num_map", NC_INT, 1, &dim_num_nodes, &var_node_num_map);
        NCERR(stat);
        stat = nc_def_var(ncid, "elem_num_map", NC_INT, 1, &dim_num_elem, &var_elem_num_map);
        NCERR(stat);

        /*----------------------------------------------------------------
         * Variables: Nodal result variable names and data  [REQUIRED]
         *----------------------------------------------------------------*/
        int var_name_nod_var = -1;
        int var_vals_nod[512]; /* enough for many variables */
        if (total_nod_var > 0) {
            int dims2[2] = {dim_num_nod_var, dim_len_name};
            stat = nc_def_var(ncid, "name_nod_var", NC_CHAR, 2, dims2, &var_name_nod_var);
            NCERR(stat);

            for (i = 0; i < total_nod_var; i++) {
                char vname[EXO_MAX_STR_LENGTH];
                snprintf(vname, EXO_MAX_STR_LENGTH, "vals_nod_var%d", i + 1);
                int dims_ts_nn[2] = {dim_time_step, dim_num_nodes};
                stat = nc_def_var(ncid, vname, NC_DOUBLE, 2, dims_ts_nn, &var_vals_nod[i]);
                NCERR(stat);
            }
        }

        /*----------------------------------------------------------------
         * Variables: Element result variable names and data  [REQUIRED]
         *----------------------------------------------------------------*/
        int var_name_elem_var = -1;
        int var_vals_elem[512]; /* [var_idx * num_blocks + blk_idx] */
        int var_elem_var_tab = -1;
        if (total_elem_var > 0) {
            int dims2[2] = {dim_num_elem_var, dim_len_name};
            stat = nc_def_var(ncid, "name_elem_var", NC_CHAR, 2, dims2, &var_name_elem_var);
            NCERR(stat);

            /* Truth table */
            {
                int dims_tt[2] = {dim_num_el_blk, dim_num_elem_var};
                stat = nc_def_var(ncid, "elem_var_tab", NC_INT, 2, dims_tt, &var_elem_var_tab);
                NCERR(stat);
            }

            /* vals_elem_var{v}eb{b} for each variable and block */
            for (i = 0; i < total_elem_var; i++) {
                for (j = 0; j < num_blocks; j++) {
                    char vname[EXO_MAX_STR_LENGTH];
                    snprintf(vname, EXO_MAX_STR_LENGTH, "vals_elem_var%deb%d", i + 1, j + 1);
                    int dims_ts_eb[2] = {dim_time_step, dim_el_in_blk[j]};
                    stat = nc_def_var(ncid, vname, NC_DOUBLE, 2, dims_ts_eb,
                                      &var_vals_elem[i * num_blocks + j]);
                    NCERR(stat);
                }
            }
        }

        /*================================================================
         * END DEFINE MODE
         *================================================================*/
        stat = nc_enddef(ncid);
        NCERR(stat);

        /*================================================================
         * WRITE STATIC DATA (coordinates, connectivity, maps, etc.)
         * These are written only once.
         *================================================================*/

        /* --- Coordinates (XYZ separated) --- */
        {
            double *cx = (double *)HECMW_malloc(sizeof(double) * n_node);
            double *cy = (double *)HECMW_malloc(sizeof(double) * n_node);
            double *cz = (double *)HECMW_malloc(sizeof(double) * n_node);
            for (i = 0; i < n_node; i++) {
                cx[i] = mesh->node[3 * i];
                cy[i] = mesh->node[3 * i + 1];
                cz[i] = mesh->node[3 * i + 2];
            }
            stat = nc_put_var_double(ncid, var_coordx, cx); NCERR(stat);
            stat = nc_put_var_double(ncid, var_coordy, cy); NCERR(stat);
            stat = nc_put_var_double(ncid, var_coordz, cz); NCERR(stat);
            HECMW_free(cx); HECMW_free(cy); HECMW_free(cz);
        }

        /* --- Coordinate names --- */
        {
            char names[3][EXO_MAX_STR_LENGTH];
            pad_string(names[0], EXO_MAX_STR_LENGTH, "x");
            pad_string(names[1], EXO_MAX_STR_LENGTH, "y");
            pad_string(names[2], EXO_MAX_STR_LENGTH, "z");
            stat = nc_put_var_text(ncid, var_coor_names, &names[0][0]);
            NCERR(stat);
        }

        /* --- Element block status and IDs --- */
        {
            int *eb_stat = (int *)HECMW_malloc(sizeof(int) * num_blocks);
            int *eb_ids  = (int *)HECMW_malloc(sizeof(int) * num_blocks);
            for (j = 0; j < num_blocks; j++) {
                eb_stat[j] = 1;      /* active */
                eb_ids[j]  = j + 1;  /* block ID = 1,2,3,... */
            }
            stat = nc_put_var_int(ncid, var_eb_status, eb_stat); NCERR(stat);
            stat = nc_put_var_int(ncid, var_eb_prop1,  eb_ids);  NCERR(stat);
            HECMW_free(eb_stat); HECMW_free(eb_ids);
        }

        /* --- [OPTIONAL] Element block names --- */
        {
            char *names_buf = (char *)HECMW_malloc(num_blocks * EXO_MAX_STR_LENGTH);
            for (j = 0; j < num_blocks; j++) {
                char label[EXO_MAX_STR_LENGTH];
                snprintf(label, EXO_MAX_STR_LENGTH, "block_%d_%s", j + 1, blocks[j].exo_name);
                pad_string(&names_buf[j * EXO_MAX_STR_LENGTH], EXO_MAX_STR_LENGTH, label);
            }
            stat = nc_put_var_text(ncid, var_eb_names, names_buf);
            NCERR(stat);
            HECMW_free(names_buf);
        }

        /* --- Connectivity per block --- */
        for (j = 0; j < num_blocks; j++) {
            int ne = blocks[j].num_elem;
            int nn = blocks[j].nod_per_elem;
            int etype = blocks[j].hecmw_type;
            int node_shift = get_node_shift(etype);
            int *conn = (int *)HECMW_malloc(sizeof(int) * ne * nn);

            for (i = 0; i < ne; i++) {
                int ei = blocks[j].elem_indices[i]; /* original element index */
                jS = mesh->elem_node_index[ei];
                jE = mesh->elem_node_index[ei + 1];

                if (etype == 342) {
                    /* TET10: reorder nodes (same table as VTK output) */
                    for (k = 0; k < nn; k++) {
                        conn[i * nn + k] = mesh->elem_node_item[jS + table342[k]];
                        /* NOTE: Exodus uses 1-based indexing (no -1) */
                    }
                } else {
                    for (k = 0; k < nn; k++) {
                        conn[i * nn + k] = mesh->elem_node_item[jS + k];
                    }
                }
            }
            stat = nc_put_var_int(ncid, var_connect[j], conn);
            NCERR(stat);
            HECMW_free(conn);
        }

        /* --- [OPTIONAL] QA records --- */
        {
            char qa_data[1][4][EXO_MAX_STR_LENGTH];
            time_t now = time(NULL);
            struct tm *t = localtime(&now);
            char datestr[32], timestr[32];
            strftime(datestr, sizeof(datestr), "%Y%m%d", t);
            strftime(timestr, sizeof(timestr), "%H:%M:%S", t);

            pad_string(qa_data[0][0], EXO_MAX_STR_LENGTH, "FrontISTR");
            pad_string(qa_data[0][1], EXO_MAX_STR_LENGTH, "1.0");
            pad_string(qa_data[0][2], EXO_MAX_STR_LENGTH, datestr);
            pad_string(qa_data[0][3], EXO_MAX_STR_LENGTH, timestr);
            stat = nc_put_var_text(ncid, var_qa_records, &qa_data[0][0][0]);
            NCERR(stat);
        }

        /* --- [OPTIONAL] Node/Element number maps --- */
        {
            int *nmap = (int *)HECMW_malloc(sizeof(int) * n_node);
            int *emap = (int *)HECMW_malloc(sizeof(int) * n_elem);
            /* Use global IDs if available, otherwise sequential */
            if (mesh->global_node_ID != NULL) {
                for (i = 0; i < n_node; i++) nmap[i] = mesh->global_node_ID[i];
            } else {
                for (i = 0; i < n_node; i++) nmap[i] = i + 1;
            }
            if (mesh->global_elem_ID != NULL) {
                for (i = 0; i < n_elem; i++) emap[i] = mesh->global_elem_ID[i];
            } else {
                for (i = 0; i < n_elem; i++) emap[i] = i + 1;
            }
            stat = nc_put_var_int(ncid, var_node_num_map, nmap); NCERR(stat);
            stat = nc_put_var_int(ncid, var_elem_num_map, emap); NCERR(stat);
            HECMW_free(nmap); HECMW_free(emap);
        }

        /* --- Nodal variable names --- */
        if (total_nod_var > 0) {
            char *name_buf = (char *)HECMW_malloc(total_nod_var * EXO_MAX_STR_LENGTH);
            int vi = 0;
            for (i = 0; i < data->nn_component; i++) {
                const char **suffixes = get_comp_suffixes(data->nn_dof[i], data->node_label[i]);
                for (k = 0; k < data->nn_dof[i]; k++) {
                    char vname[EXO_MAX_STR_LENGTH];
                    if (suffixes != NULL) {
                        snprintf(vname, EXO_MAX_STR_LENGTH, "%s%s",
                                 data->node_label[i], suffixes[k]);
                    } else {
                        snprintf(vname, EXO_MAX_STR_LENGTH, "%s_%d",
                                 data->node_label[i], k + 1);
                    }
                    pad_string(&name_buf[vi * EXO_MAX_STR_LENGTH], EXO_MAX_STR_LENGTH, vname);
                    vi++;
                }
            }
            stat = nc_put_var_text(ncid, var_name_nod_var, name_buf);
            NCERR(stat);
            HECMW_free(name_buf);
        }

        /* --- Element variable names --- */
        if (total_elem_var > 0) {
            char *name_buf = (char *)HECMW_malloc(total_elem_var * EXO_MAX_STR_LENGTH);
            int vi = 0;
            for (i = 0; i < data->ne_component; i++) {
                const char **suffixes = get_comp_suffixes(data->ne_dof[i], data->elem_label[i]);
                for (k = 0; k < data->ne_dof[i]; k++) {
                    char vname[EXO_MAX_STR_LENGTH];
                    if (suffixes != NULL) {
                        snprintf(vname, EXO_MAX_STR_LENGTH, "%s%s",
                                 data->elem_label[i], suffixes[k]);
                    } else {
                        snprintf(vname, EXO_MAX_STR_LENGTH, "%s_%d",
                                 data->elem_label[i], k + 1);
                    }
                    pad_string(&name_buf[vi * EXO_MAX_STR_LENGTH], EXO_MAX_STR_LENGTH, vname);
                    vi++;
                }
            }
            stat = nc_put_var_text(ncid, var_name_elem_var, name_buf);
            NCERR(stat);
            HECMW_free(name_buf);
        }

        /* --- Element variable truth table (all 1s: every var in every block) --- */
        if (total_elem_var > 0) {
            int *tab = (int *)HECMW_malloc(sizeof(int) * num_blocks * total_elem_var);
            for (i = 0; i < num_blocks * total_elem_var; i++) tab[i] = 1;
            stat = nc_put_var_int(ncid, var_elem_var_tab, tab);
            NCERR(stat);
            HECMW_free(tab);
        }

    } else {
        /*====================================================================
         * SUBSEQUENT CALLS (per_step==0 only): Open existing file for appending
         *====================================================================*/
        stat = nc_open(file_exo, NC_WRITE, &ncid);
        NCERR(stat);
    }

    /*========================================================================
     * WRITE TIME-DEPENDENT DATA (every call)
     *========================================================================*/

    /* --- Get time value from global data (look for TOTALTIME) --- */
    double time_value = 0.0;
    {
        int gshift = 0;
        for (i = 0; i < data->ng_component; i++) {
            if (strcmp(data->global_label[i], "TOTALTIME") == 0) {
                time_value = data->global_val_item[gshift];
                break;
            }
            gshift += data->ng_dof[i];
        }
    }

    /* Write time value */
    {
        int varid;
        /* per_step: each file has only 1 timestep (index=0)
         * normal:   timesteps accumulate in one file */
        size_t ts_idx = per_step ? 0 : (size_t)time_step_count;
        stat = nc_inq_varid(ncid, "time_whole", &varid); NCERR(stat);
        stat = nc_put_var1_double(ncid, varid, &ts_idx, &time_value); NCERR(stat);
    }

    /* --- Time step index for data arrays --- */
    size_t ts_idx = per_step ? 0 : (size_t)time_step_count;

    /* --- Write nodal variables for this time step --- */
    if (total_nod_var > 0) {
        int vi = 0;
        int comp_shift = 0;
        for (i = 0; i < data->nn_component; i++) {
            for (k = 0; k < data->nn_dof[i]; k++) {
                char vname[EXO_MAX_STR_LENGTH];
                snprintf(vname, EXO_MAX_STR_LENGTH, "vals_nod_var%d", vi + 1);
                int varid;
                stat = nc_inq_varid(ncid, vname, &varid); NCERR(stat);

                double *vals = (double *)HECMW_malloc(sizeof(double) * n_node);
                for (j = 0; j < n_node; j++) {
                    vals[j] = data->node_val_item[j * data_tot_n + comp_shift + k];
                }
                size_t start[2] = {ts_idx, 0};
                size_t count[2] = {1, (size_t)n_node};
                stat = nc_put_vara_double(ncid, varid, start, count, vals);
                NCERR(stat);
                HECMW_free(vals);
                vi++;
            }
            comp_shift += data->nn_dof[i];
        }
    }

    /* --- Write element variables for this time step --- */
    if (total_elem_var > 0) {
        int vi = 0;
        int comp_shift = 0;
        for (i = 0; i < data->ne_component; i++) {
            for (k = 0; k < data->ne_dof[i]; k++) {
                /* Write for each block */
                for (j = 0; j < num_blocks; j++) {
                    char vname[EXO_MAX_STR_LENGTH];
                    snprintf(vname, EXO_MAX_STR_LENGTH, "vals_elem_var%deb%d", vi + 1, j + 1);
                    int varid;
                    stat = nc_inq_varid(ncid, vname, &varid); NCERR(stat);

                    int ne = blocks[j].num_elem;
                    double *vals = (double *)HECMW_malloc(sizeof(double) * ne);
                    int ei_idx;
                    for (ei_idx = 0; ei_idx < ne; ei_idx++) {
                        int ei = blocks[j].elem_indices[ei_idx];
                        vals[ei_idx] = data->elem_val_item[ei * data_tot_e + comp_shift + k];
                    }
                    size_t start[2] = {ts_idx, 0};
                    size_t count[2] = {1, (size_t)ne};
                    stat = nc_put_vara_double(ncid, varid, start, count, vals);
                    NCERR(stat);
                    HECMW_free(vals);
                }
                vi++;
            }
            comp_shift += data->ne_dof[i];
        }
    }

    /*========================================================================
     * CLEANUP
     *========================================================================*/
    stat = nc_close(ncid);
    NCERR(stat);

    /* Free element block index arrays */
    for (j = 0; j < num_blocks; j++) {
        HECMW_free(blocks[j].elem_indices);
    }

    time_step_count++;
}

/*============================================================================
 * Public interface (same pattern as HECMW_vtk_output / HECMW_bin_vtk_output)
 *============================================================================*/

/* output_type=18 (EXODUS): All timesteps in one file */
void HECMW_exodus_output(struct hecmwST_local_mesh *mesh,
                         struct hecmwST_result_data *data,
                         char *outfile, char *outfile1,
                         HECMW_Comm VIS_COMM)
{
    exodus_output(mesh, data, outfile, outfile1, VIS_COMM, 0);
}

/* output_type=19 (STEP_EXODUS): One file per timestep */
void HECMW_exodus_step_output(struct hecmwST_local_mesh *mesh,
                              struct hecmwST_result_data *data,
                              char *outfile, char *outfile1,
                              HECMW_Comm VIS_COMM)
{
    exodus_output(mesh, data, outfile, outfile1, VIS_COMM, 1);
}

#else /* !WITH_NETCDF */

#include <stdio.h>
#include "hecmw_struct.h"
#include "hecmw_result.h"

/*============================================================================
 * Stub functions when compiled without NetCDF support
 *============================================================================*/
static void exodus_stub_warning(HECMW_Comm VIS_COMM)
{
    int mynode;
    HECMW_Comm_rank(VIS_COMM, &mynode);
    if (mynode == 0) {
        fprintf(stderr,
            "WARNING: output_type=EXODUS requested but FrontISTR was not\n"
            "         compiled with NetCDF support.\n"
            "         Rebuild with -DWITH_NETCDF=ON (CMake) or\n"
            "         --with-netcdf (Makefile) to enable Exodus output.\n");
    }
}

void HECMW_exodus_output(struct hecmwST_local_mesh *mesh,
                         struct hecmwST_result_data *data,
                         char *outfile, char *outfile1,
                         HECMW_Comm VIS_COMM)
{
    exodus_stub_warning(VIS_COMM);
    (void)mesh; (void)data; (void)outfile; (void)outfile1;
}

void HECMW_exodus_step_output(struct hecmwST_local_mesh *mesh,
                              struct hecmwST_result_data *data,
                              char *outfile, char *outfile1,
                              HECMW_Comm VIS_COMM)
{
    exodus_stub_warning(VIS_COMM);
    (void)mesh; (void)data; (void)outfile; (void)outfile1;
}

#endif /* WITH_NETCDF */
