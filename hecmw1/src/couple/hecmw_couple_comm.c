/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "hecmw_config.h"
#include "hecmw_msgno.h"
#include "hecmw_malloc.h"
#include "hecmw_error.h"
#include "hecmw_comm.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_info.h"
#include "hecmw_couple_comm.h"


extern int
HECMW_couple_inter_send_recv(int n_neighbor_pe_send, int *neighbor_pe_send,
		int *sendbuf_index, void *sendbuf, int n_neighbor_pe_recv, int *neighbor_pe_recv,
		int *recvbuf_index, void *recvbuf, HECMW_Datatype datatype, HECMW_Comm comm)
{
	HECMW_Request *request_send = NULL, *request_recv = NULL;
	HECMW_Status *status_send = NULL, *status_recv = NULL;
	int rtc, i;

	if(n_neighbor_pe_send > 0) {
	  request_send = (HECMW_Request *)HECMW_calloc(n_neighbor_pe_send, sizeof(HECMW_Request));
	  if(request_send == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	  status_send = (HECMW_Status *)HECMW_calloc(n_neighbor_pe_send, sizeof(HECMW_Status));
	  if(status_send == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	}
	if(n_neighbor_pe_recv > 0) {
	  request_recv = (HECMW_Request *)HECMW_calloc(n_neighbor_pe_recv, sizeof(HECMW_Request));
	  if(request_recv == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	  status_recv = (HECMW_Status *)HECMW_calloc(n_neighbor_pe_recv, sizeof(HECMW_Status));
	  if(status_recv == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	}

	if(datatype == HECMW_INT) {
		int *_sendbuf = (int *)sendbuf;
		int *_recvbuf = (int *)recvbuf;

		/* send */
		for(i=0; i<n_neighbor_pe_send; i++) {
			rtc = HECMW_Isend(&_sendbuf[sendbuf_index[i]], sendbuf_index[i+1]-sendbuf_index[i],
					HECMW_INT, neighbor_pe_send[i], 0, comm, &request_send[i]);
			if(rtc != 0)  goto error;
		}

		/* receive */
		for(i=0; i<n_neighbor_pe_recv; i++) {
			rtc = HECMW_Irecv(&_recvbuf[recvbuf_index[i]], recvbuf_index[i+1]-recvbuf_index[i],
					HECMW_INT, neighbor_pe_recv[i], 0, comm, &request_recv[i]);
			if(rtc != 0)  goto error;
		}
	} else if(datatype == HECMW_DOUBLE) {
		double *_sendbuf = (double *)sendbuf;
		double *_recvbuf = (double *)recvbuf;

		/* send */
		for(i=0; i<n_neighbor_pe_send; i++) {
			rtc = HECMW_Isend(&_sendbuf[sendbuf_index[i]], sendbuf_index[i+1]-sendbuf_index[i],
					HECMW_DOUBLE, neighbor_pe_send[i], 0, comm, &request_send[i]);
			if(rtc != 0)  goto error;
		}

		/* receive */
		for(i=0; i<n_neighbor_pe_recv; i++) {
			rtc = HECMW_Irecv(&_recvbuf[recvbuf_index[i]], recvbuf_index[i+1]-recvbuf_index[i],
					HECMW_DOUBLE, neighbor_pe_recv[i], 0, comm, &request_recv[i]);
			if(rtc != 0)  goto error;
		}
	} else if(datatype == HECMW_CHAR) {
		char *_sendbuf = (char *)sendbuf;
		char *_recvbuf = (char *)recvbuf;

		/* send */
		for(i=0; i<n_neighbor_pe_send; i++) {
			rtc = HECMW_Isend(&_sendbuf[sendbuf_index[i]], sendbuf_index[i+1]-sendbuf_index[i],
					HECMW_CHAR, neighbor_pe_send[i], 0, comm, &request_send[i]);
			if(rtc != 0)  goto error;
		}

		/* receive */
		for(i=0; i<n_neighbor_pe_recv; i++) {
			rtc = HECMW_Irecv(&_recvbuf[recvbuf_index[i]], recvbuf_index[i+1]-recvbuf_index[i],
					HECMW_CHAR, neighbor_pe_recv[i], 0, comm, &request_recv[i]);
			if(rtc != 0)  goto error;
		}
	} else {
		HECMW_set_error(HECMWCPL_E_MPI_DATATYPE, "");
		goto error;
	}

	/* wait */
	if(n_neighbor_pe_recv > 0) {
		rtc = HECMW_Waitall(n_neighbor_pe_recv, request_recv, status_recv);
		if(rtc != 0)  goto error;
	}

	if(n_neighbor_pe_send > 0) {
		rtc = HECMW_Waitall(n_neighbor_pe_send, request_send, status_send);
		if(rtc != 0)  goto error;
	}

	HECMW_free(request_send);
	HECMW_free(request_recv);
	HECMW_free(status_send);
	HECMW_free(status_recv);

	return 0;

error:
	HECMW_free(request_send);
	HECMW_free(request_recv);
	HECMW_free(status_send);
	HECMW_free(status_recv);

	return -1;
}


extern int
HECMW_couple_intra_send_recv(int n_neighbor_pe, int *neighbor_pe, int *sendbuf_index,
		void *sendbuf, int *recvbuf_index, void *recvbuf, HECMW_Datatype datatype, HECMW_Comm comm)
{
	HECMW_Request *request_send = NULL, *request_recv = NULL;
	HECMW_Status *status_send = NULL, *status_recv = NULL;
	int rtc, i;

	if(n_neighbor_pe > 0) {
	  request_send = (HECMW_Request *)HECMW_calloc(n_neighbor_pe, sizeof(HECMW_Request));
	  if(request_send == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	  status_send = (HECMW_Status *)HECMW_calloc(n_neighbor_pe, sizeof(HECMW_Status));
	  if(status_send == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	  request_recv = (HECMW_Request *)HECMW_calloc(n_neighbor_pe, sizeof(HECMW_Request));
	  if(request_recv == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	  status_recv = (HECMW_Status *)HECMW_calloc(n_neighbor_pe, sizeof(HECMW_Status));
	  if(status_recv == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	}

	if(datatype == HECMW_INT) {
		int *_sendbuf = (int *)sendbuf;
		int *_recvbuf = (int *)recvbuf;

		/* send */
		for(i=0; i<n_neighbor_pe; i++) {
			rtc = HECMW_Isend(&_sendbuf[sendbuf_index[i]], sendbuf_index[i+1]-sendbuf_index[i],
					datatype, neighbor_pe[i], 0, comm, &request_send[i]);
			if(rtc != 0)  goto error;
		}

		/* receive */
		for(i=0; i<n_neighbor_pe; i++) {
			rtc = HECMW_Irecv(&_recvbuf[recvbuf_index[i]], recvbuf_index[i+1]-recvbuf_index[i],
					datatype, neighbor_pe[i], 0, comm, &request_recv[i]);
			if(rtc != 0)  goto error;
		}
	} else if(datatype == HECMW_DOUBLE) {
		double *_sendbuf = (double *)sendbuf;
		double *_recvbuf = (double *)recvbuf;

		/* send */
		for(i=0; i<n_neighbor_pe; i++) {
			rtc = HECMW_Isend(&_sendbuf[sendbuf_index[i]], sendbuf_index[i+1]-sendbuf_index[i],
					datatype, neighbor_pe[i], 0, comm, &request_send[i]);
			if(rtc != 0)  goto error;
		}

		/* receive */
		for(i=0; i<n_neighbor_pe; i++) {
			rtc = HECMW_Irecv(&_recvbuf[recvbuf_index[i]], recvbuf_index[i+1]-recvbuf_index[i],
					datatype, neighbor_pe[i], 0, comm, &request_recv[i]);
			if(rtc != 0)  goto error;
		}
	} else if(datatype == HECMW_CHAR) {
		char *_sendbuf = (char *)sendbuf;
		char *_recvbuf = (char *)recvbuf;

		/* send */
		for(i=0; i<n_neighbor_pe; i++) {
			rtc = HECMW_Isend(&_sendbuf[sendbuf_index[i]], sendbuf_index[i+1]-sendbuf_index[i],
					datatype, neighbor_pe[i], 0, comm, &request_send[i]);
			if(rtc != 0)  goto error;
		}

		/* receive */
		for(i=0; i<n_neighbor_pe; i++) {
			rtc = HECMW_Irecv(&_recvbuf[recvbuf_index[i]], recvbuf_index[i+1]-recvbuf_index[i],
					datatype, neighbor_pe[i], 0, comm, &request_recv[i]);
			if(rtc != 0)  goto error;
		}
	} else {
		HECMW_set_error(HECMWCPL_E_MPI_DATATYPE, "");
		goto error;
	}

	/* wait */
	if(n_neighbor_pe > 0) {
		rtc = HECMW_Waitall(n_neighbor_pe, request_recv, status_recv);
		if(rtc != 0)  goto error;
	}

	if(n_neighbor_pe > 0) {
		rtc = HECMW_Waitall(n_neighbor_pe, request_send, status_send);
		if(rtc != 0)  goto error;
	}

	HECMW_free(request_send);
	HECMW_free(request_recv);
	HECMW_free(status_send);
	HECMW_free(status_recv);

	return 0;

error:
	HECMW_free(request_send);
	HECMW_free(request_recv);
	HECMW_free(status_send);
	HECMW_free(status_recv);

	return -1;
}


extern int
HECMW_couple_bcast(int n_neighbor_pe_send, int *neighbor_pe_send,
		int sendbuf_size, void *sendbuf, int n_neighbor_pe_recv, int *neighbor_pe_recv,
		int *recvbuf_index, void *recvbuf, HECMW_Datatype datatype, HECMW_Comm comm)
{
	HECMW_Request *request_send = NULL, *request_recv = NULL;
	HECMW_Status *status_send = NULL, *status_recv = NULL;
	int rtc, i;

	if(n_neighbor_pe_send > 0) {
	  request_send = (HECMW_Request *)HECMW_calloc(n_neighbor_pe_send, sizeof(HECMW_Request));
	  if(request_send == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	  status_send = (HECMW_Status *)HECMW_calloc(n_neighbor_pe_send, sizeof(HECMW_Status));
	  if(status_send == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	}
	if(n_neighbor_pe_recv > 0) {
	  request_recv = (HECMW_Request *)HECMW_calloc(n_neighbor_pe_recv, sizeof(HECMW_Request));
	  if(request_recv == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	  status_recv = (HECMW_Status *)HECMW_calloc(n_neighbor_pe_recv, sizeof(HECMW_Status));
	  if(status_recv == NULL) {
		  HECMW_set_error(errno, "");
		  goto error;
	  }
	}

	if(datatype == HECMW_INT) {
		int *_sendbuf = (int *)sendbuf;
		int *_recvbuf = (int *)recvbuf;

		/* send */
		for(i=0; i<n_neighbor_pe_send; i++) {
			rtc = HECMW_Isend(&_sendbuf[0], sendbuf_size,
					datatype, neighbor_pe_send[i], 0, comm, &request_send[i]);
			if(rtc != 0)  goto error;
		}

		/* receive */
		for(i=0; i<n_neighbor_pe_recv; i++) {
			rtc = HECMW_Irecv(&_recvbuf[recvbuf_index[i]], recvbuf_index[i+1]-recvbuf_index[i],
					datatype, neighbor_pe_recv[i], 0, comm, &request_recv[i]);
			if(rtc != 0)  goto error;
		}
	} else if(datatype == HECMW_DOUBLE) {
		double *_sendbuf = (double *)sendbuf;
		double *_recvbuf = (double *)recvbuf;

		/* send */
		for(i=0; i<n_neighbor_pe_send; i++) {
			rtc = HECMW_Isend(&_sendbuf[0], sendbuf_size,
					datatype, neighbor_pe_send[i], 0, comm, &request_send[i]);
			if(rtc != 0)  goto error;
		}

		/* receive */
		for(i=0; i<n_neighbor_pe_recv; i++) {
			rtc = HECMW_Irecv(&_recvbuf[recvbuf_index[i]], recvbuf_index[i+1]-recvbuf_index[i],
					datatype, neighbor_pe_recv[i], 0, comm, &request_recv[i]);
			if(rtc != 0)  goto error;
		}
	} else if(datatype == HECMW_CHAR) {
		char *_sendbuf = (char *)sendbuf;
		char *_recvbuf = (char *)recvbuf;

		/* send */
		for(i=0; i<n_neighbor_pe_send; i++) {
			rtc = HECMW_Isend(&_sendbuf[0], sendbuf_size,
					datatype, neighbor_pe_send[i], 0, comm, &request_send[i]);
			if(rtc != 0)  goto error;
		}

		/* receive */
		for(i=0; i<n_neighbor_pe_recv; i++) {
			rtc = HECMW_Irecv(&_recvbuf[recvbuf_index[i]], recvbuf_index[i+1]-recvbuf_index[i],
					datatype, neighbor_pe_recv[i], 0, comm, &request_recv[i]);
			if(rtc != 0)  goto error;
		}
	} else {
		HECMW_set_error(HECMWCPL_E_MPI_DATATYPE, "");
		goto error;
	}

	/* wait */
	if(n_neighbor_pe_recv > 0) {
		rtc = HECMW_Waitall(n_neighbor_pe_recv, request_recv, status_recv);
		if(rtc != 0)  goto error;
	}
	if(n_neighbor_pe_send > 0) {
		rtc = HECMW_Waitall(n_neighbor_pe_send, request_send, status_send);
		if(rtc != 0)  goto error;
	}

	HECMW_free(request_send);
	HECMW_free(request_recv);
	HECMW_free(status_send);
	HECMW_free(status_recv);

	return 0;

error:
	HECMW_free(request_send);
	HECMW_free(request_recv);
	HECMW_free(status_send);
	HECMW_free(status_recv);
	return -1;
}
#if 0
/*================================================================================================*/

static int
send_recv_r2r_int(const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_comm *comm_dst, int *buffer, int count)
{
	int n_pe_send = 0, n_pe_recv = 0, *pe_send = NULL, *pe_recv = NULL;
	int *sendbuf_index = NULL, *recvbuf_index = NULL, *sendbuf = NULL, *recvbuf = NULL;
	int size, rtc, i;

	if(comm_src->is_root) {
		n_pe_send = 1;

		pe_send = (int *)HECMW_malloc(sizeof(int)*n_pe_send);
		if(pe_send == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		pe_send[0] = comm_dst->root;

		sendbuf_index = (int *)HECMW_calloc(n_pe_send+1, sizeof(int));
		if(sendbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		sendbuf_index[1] = count;

		sendbuf = (int *)HECMW_malloc(sizeof(int)*sendbuf_index[1]);
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<sendbuf_index[1]; i++) {
			sendbuf[i] = buffer[i];
		}
	}

	if(comm_dst->is_root) {
		n_pe_recv = 1;

		pe_recv = (int *)HECMW_malloc(sizeof(int)*n_pe_recv);
		if(pe_recv == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		pe_recv[0] = comm_src->root;

		recvbuf_index = (int *)HECMW_calloc(n_pe_recv+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		recvbuf_index[1] = count;

		recvbuf = (int *)HECMW_malloc(sizeof(int)*recvbuf_index[1]);
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	rtc = HECMW_couple_inter_send_recv(n_pe_send, pe_send, sendbuf_index, sendbuf,
			n_pe_recv, pe_recv, recvbuf_index, recvbuf, HECMW_INT, HECMW_comm_get_comm());
	if(rtc) goto error;

	if(comm_dst->is_root) {
		memcpy(buffer, recvbuf, size);
	}
return 0;
	HECMW_free(pe_send);
	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(pe_recv);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	return 0;

error:
	HECMW_free(pe_send);
	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(pe_recv);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	return -1;
}
#endif
#if 0

extern int
HECMW_couple_inter_bcast(char *boundary_id, void *sendbuf, int sendcount,
		void **recvbuf, int recvcount, HECMW_Datatype datatype, int direction)
{
	struct hecmw_couple_comm *comm_src = NULL, *comm_dst = NULL, *intercomm = NULL;
	int rtc, i;

	if((intercomm = HECMW_couple_get_intercomm(boundary_id)) == NULL) goto error;

	if(direction == HECMW_COUPLE_UNIT1_TO_UNIT2) {
		comm_src = HECMW_couple_get_intracomm(boundary_id, HECMW_COUPLE_UNIT1);
		comm_dst = HECMW_couple_get_intracomm(boundary_id, HECMW_COUPLE_UNIT2);
	} else if(direction == HECMW_COUPLE_UNIT2_TO_UNIT1) {
		comm_src = HECMW_couple_get_intracomm(boundary_id, HECMW_COUPLE_UNIT2);
		comm_dst = HECMW_couple_get_intracomm(boundary_id, HECMW_COUPLE_UNIT1);
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_DIRECTION, "");
		goto error;
	}
	if(comm_src == NULL || comm_dst == NULL) goto error;

	if(datatype == HECMW_INT) {
		int n_src = 0, n_dst = 0, *src = NULL, *dst = NULL;
		int *_sendbuf = (int *)sendbuf;
		int *sendindex = NULL, *recvindex = NULL;
		int _recvcount, *_recvbuf = NULL;

		if(comm_src->is_root) {
			n_dst = 1;

			dst = (int *)HECMW_malloc(sizeof(int)*n_dst);
			if(dst == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			dst[0] = comm_dst->root;

			sendindex = (int *)HECMW_calloc(n_dst+1, sizeof(int));
			if(sendindex == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			sendindex[0] = 0;
			sendindex[1] = 1;

			_sendbuf = (int *)HECMW_malloc(sizeof(int)*sendindex[n_dst]);
			if(_sendbuf == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			_sendbuf[0] = sendcount;
		}

		if(comm_dst->is_root) {
			n_src = 1;

			src = (int *)HECMW_malloc(sizeof(int)*n_src);
			if(src == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			src[0] = comm_src->root;

			recvindex = (int *)HECMW_calloc(n_src+1, sizeof(int));
			if(recvindex == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			recvindex[0] = 0;
			recvindex[1] = 1;

			_recvbuf = (int *)HECMW_malloc(sizeof(int)*recvindex[n_src]);
			if(_recvbuf == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
		}

		rtc = HECMW_couple_inter_send_recv(n_dst, dst, sendindex, _sendbuf,
				n_src, src, recvindex, _recvbuf, HECMW_INT, intercomm->comm);
		if(rtc) goto error;

		if(comm_dst->is_root) {
			recvcount = _recvbuf[0];
		}

		HECMW_free(_sendbuf);
		HECMW_free(sendindex);
		HECMW_free(_recvbuf);
		HECMW_free(recvindex);
		_sendbuf = NULL;
		sendindex = NULL;
		_recvbuf = NULL;
		recvindex = NULL;


		if(comm_src->is_root) {
			n_dst = 1;

			dst = (int *)HECMW_malloc(sizeof(int)*n_dst);
			if(dst == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			dst[0] = comm_dst->root;

			sendindex = (int *)HECMW_calloc(n_dst+1, sizeof(int));
			if(sendindex == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			sendindex[1] = sendcount;

			_sendbuf = (int *)HECMW_malloc(sizeof(int)*sendindex[n_dst]+1);
			if(_sendbuf == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(i=0; i<sendindex[n_dst]; i++) {
				_sendbuf[i] = sendbuf[i];
			}
		}

		if(comm_dst->is_root) {
			n_src = 1;

			src = (int *)HECMW_malloc(sizeof(int)*n_src);
			if(src == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			src[0] = comm_src->root;

			recvindex = (int *)HECMW_calloc(n_src+1, sizeof(int));
			if(recvindex == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			recvindex[1] = recvcount;

			_recvbuf = (int *)HECMW_malloc(sizeof(int)*recvindex[n_src]);
			if(_recvbuf == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
		}

		rtc = HECMW_couple_inter_send_recv(n_dst, dst, sendindex, _sendbuf,
				n_src, src, recvindex, _recvbuf, HECMW_INT, intercomm->comm);
		if(rtc) goto error;

		if(comm_dst->is_member) {
			if(HECMW_Bcast(recvcount, 1, HECMW_INT, 0,  comm_dst->comm)) goto error;

			recvbuf = HECMW_malloc(sizeof(int)*recvcount);
			if(recvbuf == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}

			if(comm_dst->is_root) {
				for(i=0; i<recvcount; i++) {
					recvbuf[i] = _recvbuf[i];
				}
			}
			if(HECMW_Bcast(recvbuf, recvcount, HECMW_INT, 0, comm_dst->comm)) goto error;
		}
	}

	return 0;

error:
	HECMW_couple_free_comm(intercomm);
	HECMW_couple_free_comm(comm_src);
	HECMW_couple_free_comm(comm_dst);
	return -1;
}



extern int
HECMW_couple_inter_barrier(const char *boundary_id)
{
	struct hecmw_couple_comm *intercomm = NULL;
	int rtc;

	if((intercomm = HECMW_couple_get_intercomm(boundary_id)) == NULL) return -1;

	rtc = HECMW_Barrier(intercomm->comm);
	if(rtc != HECMW_SUCCESS) goto error;

	HECMW_couple_free_comm(intercomm);

	return 0;

error:
	HECMW_couple_free_comm(intercomm);
	return -1;
}
#endif
