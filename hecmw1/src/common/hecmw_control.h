/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_CONTROL_INCLUDED
#define HECMW_CONTROL_INCLUDED

#define HECMW_CTRL_FILE "hecmw_ctrl.dat"

struct hecmw_ctrl_meshfile {
  int type;

#define HECMW_CTRL_FTYPE_HECMW_DIST 1

#define HECMW_CTRL_FTYPE_HECMW_ENTIRE 2

#define HECMW_CTRL_FTYPE_GEOFEM 3

#define HECMW_CTRL_FTYPE_ABAQUS 4

#define HECMW_CTRL_FTYPE_NASTRAN 5

#define HECMW_CTRL_FTYPE_FEMAP 6

  int io;

#define HECMW_CTRL_FILE_IO_IN 1

#define HECMW_CTRL_FILE_IO_OUT 2

  int refine;

  char *filename;
};

#define HECMW_CTRL_FILE_IO_INOUT 4

struct hecmw_ctrl_meshfiles {
  int n_mesh;

  struct hecmw_ctrl_meshfile *meshfiles;
};

extern int HECMW_ctrl_init(void);
extern int HECMW_ctrl_init_ex(const char *ctrlfile);
extern int HECMW_ctrl_finalize(void);
extern struct hecmw_ctrl_meshfiles *HECMW_ctrl_get_meshfiles(char *name_ID);
extern struct hecmw_ctrl_meshfiles *HECMW_ctrl_get_meshfiles_header(
    char *name_ID);
extern struct hecmw_ctrl_meshfiles *HECMW_ctrl_get_meshfiles_sub(char *name_ID,
                                                                 int n_rank,
                                                                 int i_rank);
extern struct hecmw_ctrl_meshfiles *HECMW_ctrl_get_meshfiles_header_sub(
    char *name_ID, int n_rank, int i_rank);
extern void HECMW_ctrl_free_meshfiles(struct hecmw_ctrl_meshfiles *meshfiles);
extern char *HECMW_ctrl_get_result_file(char *name_ID, int istep,
                                        int *fg_text);
extern char *HECMW_ctrl_get_result_fileheader(char *name_ID,
                                              int istep, int *fg_text);
extern char *HECMW_ctrl_get_result_file_sub(char *name_ID, int istep,
                                            int n_rank, int i_rank,
                                            int *fg_text);
extern char *HECMW_ctrl_get_result_fileheader_sub(char *name_ID,
                                                  int istep, int n_rank,
                                                  int i_rank, int *fg_text);
extern char *HECMW_ctrl_get_result_filebody(char *name_ID);
extern char *HECMW_ctrl_get_restart_file(char *name_ID);
extern char *HECMW_ctrl_get_restart_file_by_io(int io);
extern char *HECMW_ctrl_get_control_file(char *name_ID);
extern int HECMW_ctrl_is_exists_control(char *name_ID);
extern int HECMW_ctrl_make_subdir(char *filename);
extern int HECMW_ctrl_is_subdir(void);

#endif
