/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_COMM
#define INC_HECMW_COUPLE_COMM

#include "hecmw_config.h"

extern int HECMW_couple_inter_send_recv(
    int n_neighbor_pe_send, int *neighbor_pe_send, int *sendbuf_index,
    void *sendbuf, int n_neighbor_pe_recv, int *neighbor_pe_recv,
    int *recvbuf_index, void *recvbuf, HECMW_Datatype datatype,
    HECMW_Comm comm);

extern int HECMW_couple_intra_send_recv(int n_neighbor_pe, int *neighbor_pe,
                                        int *sendbuf_index, void *sendbuf,
                                        int *recvbuf_index, void *recvbuf,
                                        HECMW_Datatype datatype,
                                        HECMW_Comm comm);

extern int HECMW_couple_bcast(int n_neighbor_pe_send, int *neighbor_pe_send,
                              int sendbuf_size, void *sendbuf,
                              int n_neighbor_pe_recv, int *neighbor_pe_recv,
                              int *recvbuf_index, void *recvbuf,
                              HECMW_Datatype datatype, HECMW_Comm comm);

extern int HECMW_couple_inter_bcast(const char *boundary_id, void *buffer,
                                    int count, HECMW_Datatype datatype,
                                    int direction);

extern int HECMW_couple_inter_barrier(const char *boundary_id);

#endif /* INC_HECMW_COUPLE_COMM */
