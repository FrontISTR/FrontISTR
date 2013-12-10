/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Dynamic Load Balancing                            *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

#include "hecmw_repart.h"

/********** MOVE TO hecmw_comm.c (05/08/10 by N.Imai) ***********

int
HECMW_Comm_rank( HECMW_Comm comm, int *rank )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Comm_rank( comm, rank );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Comm_rank" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	    *rank=0;
        return 0;
#endif
}

int
HECMW_Comm_size( HECMW_Comm comm, int *size )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Comm_size( comm, size );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Comm_size" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	    *size=1;
        return 0;
#endif
}

int
HECMW_Comm_dup( HECMW_Comm comm, HECMW_Comm *new_comm )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Comm_dup(comm, new_comm);
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Comm_dup" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	    *new_comm=0;
        return 0;
#endif
}


int
HECMW_Finalize()
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Finalize();
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Finalize" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	    exit(0);
#endif
}

***********************************************************/


void whole_copy_array(int *recv_num, int *global_recv_num, int mynode, int pesize, HECMW_Comm repart_comm)
{
	int  i,j;
	int *tmp_recv;
	HECMW_Status  stat;

	if(mynode==0) {
		for(j=0;j<pesize+1;j++)
			global_recv_num[j]=recv_num[j];
		tmp_recv=(int *)calloc(pesize+1, sizeof(int));
		if(tmp_recv==NULL)
			HECMW_vis_memory_exit("tmp_recv");
		for(i=1;i<pesize;i++) {
			HECMW_Recv(tmp_recv, pesize+1, HECMW_INT, i, HECMW_ANY_TAG, repart_comm, &stat);
			for(j=0;j<pesize+1;j++)
				global_recv_num[i*(pesize+1)+j]=tmp_recv[j];
		}
		for(i=1;i<pesize;i++)
	       HECMW_Send(global_recv_num,(pesize+1)*pesize,HECMW_INT, i, 0, repart_comm);
		free(tmp_recv);
	}
	else {
		HECMW_Send(recv_num, pesize+1, HECMW_INT, 0, 0, repart_comm);
        HECMW_Recv(global_recv_num, (pesize+1)*pesize, HECMW_INT, 0, HECMW_ANY_TAG, repart_comm, &stat);
	}
	return;
}





int  stack_part_send_recv(int neibpetot, int *neibpe, int *stack_import,  int *stack_export, 
		     HECMW_Comm repart_comm, int my_rank)                                      
{

	  HECMW_Status	*sta1, *sta2;
  HECMW_Request	*req1, *req2;

  int	nflag = 0;
  int	neib;
  int	istart,  num;
  int	k;

  if (nflag == 0) {
    sta1 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status));
    if (sta1 == NULL) 
		HECMW_vis_memory_exit("HECMW_STATUS: stat1");
    sta2 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status));
    if (sta2 == NULL) 
		HECMW_vis_memory_exit("HECMW_STATUS: stat2");
    if ((req1 = (HECMW_Request *)calloc(neibpetot, sizeof(HECMW_Request)))
	== NULL) 
	    HECMW_vis_memory_exit("HECMW_STATUS: req1");
    if ((req2 = (HECMW_Request *)calloc(neibpetot, sizeof(HECMW_Request)))
	== NULL) 
	   HECMW_vis_memory_exit("HECMW_STATUS: req2");
    nflag = 1;
  }

  for(neib=0;neib<neibpetot;neib++) {
	  num=stack_import[neib];
    HECMW_Isend(&num, 1, HECMW_INT, neibpe[neib], 0, repart_comm, &req1[neib]);
  }
  for (neib = 0; neib < neibpetot; neib++) {
    HECMW_Irecv(&stack_export[neib], 1, HECMW_INT,neibpe[neib], 0, repart_comm, &req2[neib]);
/*	fprintf(stderr, "PE %d recv %d from %d\n", my_rank, stack_export[neib], neib);
*/
  }
  HECMW_Barrier(repart_comm);
  free(sta1);
  free(sta2);
  free(req1);
  free(req2);
       
   return;
}

int  stack_whole_send_recv(int pesize, int *stack_export,int *stack_import, HECMW_Comm repart_comm, int my_rank)                                    
{
  HECMW_Status stat;
  int tmp_int;
  int i,j;

/*
	  HECMW_Status	*sta1, *sta2;
  HECMW_Request	*req1, *req2;

  int	nflag = 0;
  int	neib;
  int	istart, num, num1;
  int	k;

  stack_export[0]=0;
  stack_export[my_rank+1]=stack_import[my_rank+1]-stack_import[my_rank];

  if (nflag == 0) {
    sta1 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE*(pesize-1), sizeof(HECMW_Status));
    if (sta1 == NULL) 
		HECMW_vis_memory_exit("HECMW_STATUS: stat1");
    sta2 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE*(pesize-1), sizeof(HECMW_Status));
    if (sta2 == NULL) 
		HECMW_vis_memory_exit("HECMW_STATUS: stat2");
    if ((req1 = (HECMW_Request *)calloc(pesize-1, sizeof(HECMW_Request)))
	== NULL) 
	    HECMW_vis_memory_exit("HECMW_STATUS: req1");
    if ((req2 = (HECMW_Request *)calloc(pesize-1, sizeof(HECMW_Request)))
	== NULL) 
	   HECMW_vis_memory_exit("HECMW_STATUS: req2");
    nflag = 1;
  }

  for(neib=0;neib<pesize;neib++) {
	  if(neib!=my_rank) {
	  num=stack_import[neib+1]-stack_import[neib];
    HECMW_Isend(&num, 1, HECMW_INT,
	      neib, 0, repart_comm, &req1[neib]);
	fprintf(stderr, "pe %d send %d to %d\n", my_rank, num, neib);
  }
  }
  for (neib = 0; neib < pesize; neib++) {
	  if(neib!=my_rank) {
    HECMW_Irecv(&stack_export[neib+1], 1, HECMW_INT,
	      neib, 0, repart_comm, &req2[neib]);

	fprintf(stderr, "pe %d recv %d from %d\n", my_rank, stack_export[neib+1], neib);


	  }
  }
		  

  HECMW_Barrier(repart_comm);

  free(sta1);
  free(sta2);
  free(req1);
  free(req2);
  */
		for(i=0;i<pesize;i++) {
		if(i!=my_rank) {
			tmp_int=stack_export[i+1]-stack_export[i];
		   HECMW_Send(&tmp_int, 1,HECMW_INT, i, 0, repart_comm);
		}
		else if(i==my_rank) {
			stack_import[i+1]=stack_export[i+1]-stack_export[i];
			for(j=0;j<pesize;j++) {
				if(j!=my_rank) {
					  HECMW_Recv(&stack_import[j+1], 1, HECMW_INT, j, HECMW_ANY_TAG, repart_comm, &stat);
				}
			}
			stack_import[0]=0;
			for(j=1;j<pesize+1;j++) 
				stack_import[j]=stack_import[j]+stack_import[j-1];
		}
	}

   HECMW_Barrier(repart_comm);
       
   return;
}


int int_part_send_recv(int n, int neibpetot, int *neibpe,int *stack_import, int *nod_import,int *stack_export, int *nod_export, 
		     int *x,  HECMW_Comm repart_comm, int my_rank)
{
  /* Important:: node ID in nod_import and nod_export are all starting from 1  */

  HECMW_Status	*sta1, *sta2;
  HECMW_Request	*req1, *req2;

  int	nflag = 0;
  int	neib;
  int	istart, inum;
  int	k;
  int   *ws, *wr;

  ws=(int *)calloc(n, sizeof(int));
  wr=(int *)calloc(n, sizeof(int));
  if((ws==NULL) || (wr==NULL)) 
	  HECMW_vis_memory_exit("send_recv: ws, wr");
  if (nflag == 0) {
    sta1 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status));
    if (sta1 == NULL) 
		HECMW_vis_memory_exit("HECMW_STATUS: stat1");
    sta2 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status));
    if (sta2 == NULL) 
		HECMW_vis_memory_exit("HECMW_STATUS: stat2");
    if ((req1 = (HECMW_Request *)calloc(neibpetot, sizeof(HECMW_Request)))
	== NULL) 
	    HECMW_vis_memory_exit("HECMW_STATUS: req1");
    if ((req2 = (HECMW_Request *)calloc(neibpetot, sizeof(HECMW_Request)))
	== NULL) 
	   HECMW_vis_memory_exit("HECMW_STATUS: req2");
    nflag = 1;
  }
  /* SEND */
  for (neib = 0; neib < neibpetot; neib++) {
/*    if (neib != 0) istart = stack_export[neib - 1];
    else istart = 0;
*/
    inum   = stack_export[neib+1] - stack_export[neib];
    
    for (k = stack_export[neib]; k < stack_export[neib] + inum; k++) {
      ws[k] = x[nod_export[k]-1];
    }
    HECMW_Isend(&ws[stack_export[neib]], inum, HECMW_INT,
	      neibpe[neib], 0, repart_comm, &req1[neib]);
  }

  /* RECEIVE */
  for (neib = 0; neib < neibpetot; neib++) {

   inum = stack_import[neib+1] - stack_import[neib];
    HECMW_Irecv(&wr[stack_import[neib]], inum, HECMW_INT,
	      neibpe[neib], 0, repart_comm, &req2[neib]);
  }
  HECMW_Barrier(repart_comm);
  /*
  HECMW_Waitall(neibpetot, req2, sta2);
  */
  for (neib = 0; neib < neibpetot; neib++) {
/*    if (neib != 0) istart = stack_import[neib-1];
    else istart = 0;
*/
    inum = stack_import[neib+1] - stack_import[neib];
    for (k = stack_import[neib]; k < stack_import[neib]+inum; k++) {
      x[nod_import[k]-1] = wr[k];
    }
  }
  HECMW_Barrier(repart_comm);
  /*
  HECMW_Waitall(neibpetot, req1, sta1);
  */
  free(sta1);
  free(sta2);
  free(req1);
  free(req2);
  free(ws);
  free(wr);
  return 1;
}


int double_part_send_recv(int n, int neibpetot, int *neibpe, int *stack_import, int *nod_import,
		     int *stack_export, int *nod_export, double *x, HECMW_Comm repart_comm, int my_rank)
{
  HECMW_Status	*sta1, *sta2;
  HECMW_Request	*req1, *req2;

  int	nflag = 0;
  int	neib;
  int	istart, inum;
  int	k;
  double   *ws, *wr;

  ws=(double *)calloc(n, sizeof(double));
  wr=(double *)calloc(n, sizeof(double));
  if((ws==NULL) || (wr==NULL)) 
	  HECMW_vis_memory_exit("send_recv: ws, wr");
  if (nflag == 0) {
    sta1 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE*neibpetot, sizeof(HECMW_Status));
    if (sta1 == NULL) 
		HECMW_vis_memory_exit("send_recv: sta1");
    sta2 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE*neibpetot, sizeof(HECMW_Status));
    if (sta2 == NULL) 
		HECMW_vis_memory_exit("send_recv: sta12");
    if ((req1 = (HECMW_Request *)calloc(neibpetot, sizeof(HECMW_Request)))
	== NULL) 
	 HECMW_vis_memory_exit("send_recv: req1");
    if ((req2 = (HECMW_Request *)calloc(neibpetot, sizeof(HECMW_Request)))
	== NULL) 
	HECMW_vis_memory_exit("send_recv: req2");
    nflag = 1;
  }
  /* SEND */
  for (neib = 0; neib < neibpetot; neib++) {
/*    if (neib != 0) istart = stack_export[neib - 1];
    else istart = 0;
*/
    inum   = stack_export[neib+1] - stack_export[neib];
    
    for (k = stack_export[neib]; k < stack_export[neib] + inum; k++) {
      ws[k] = x[nod_export[k]-1];
    }
	if(inum>0) {
/*	fprintf(stderr, "PE %d send %d node to PE %d\n", my_rank, inum, neibpe[neib]);
*/
    HECMW_Isend(&ws[stack_export[neib]], inum, HECMW_DOUBLE,
	      neibpe[neib], 0, repart_comm, &req1[neib]);
	}
  }

  /* RECEIVE */
  for (neib = 0; neib < neibpetot; neib++) {

   inum = stack_import[neib+1] - stack_import[neib];
   if(inum>0) {
/*		fprintf(stderr, "PE %d recieve %d node from PE %d\n", my_rank, inum, neibpe[neib]);
*/
    HECMW_Irecv(&wr[stack_import[neib]], inum, HECMW_DOUBLE,
	      neibpe[neib], 0, repart_comm, &req2[neib]);
   }
  }
  HECMW_Barrier(repart_comm);
/*  HECMW_Waitall(neibpetot, req2, sta2);
*/
  for (neib = 0; neib < neibpetot; neib++) {
/*    if (neib != 0) istart = stack_import[neib-1];
    else istart = 0;
*/
    inum = stack_import[neib+1] - stack_import[neib];
    for (k = stack_import[neib]; k < stack_import[neib]+inum; k++) {
      x[nod_import[k]-1] = wr[k];
    }
  }
   HECMW_Barrier(repart_comm);

/*  HECMW_Waitall(neibpetot, req1, sta1);
  */
  free(sta1);
  free(sta2);
  free(req1);
  free(req2);
  free(ws);
  free(wr);
  return 1;
}







void int_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, int *nod_import,
		     int *stack_export, int *nod_export, 
		     int *x, int *y, 
		     HECMW_Comm repart_comm, int my_rank)
{
  /* Important:: node ID in nod_import and nod_export are all starting from 0  */

/*  HECMW_Status	*sta1, *sta2;
*/
  HECMW_Request	*req1, *req2;

  int	nflag = 0;
  int	neib;
  int	istart, inum;
  int	k;
  int   *ws, *wr;

  ws=(int *)calloc(n1, sizeof(int));
  wr=(int *)calloc(n2, sizeof(int));
  if((ws==NULL) || (wr==NULL)) 
	  HECMW_vis_memory_exit("ws, wr");
  if (nflag == 0) {
/*    sta1 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status));
    if (sta1 == NULL) {
      fprintf(stderr, "Not enough memory\n");
      exit(1);
    }
    sta2 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status));
    if (sta2 == NULL) {
      fprintf(stderr, "Not enough memory\n");
      exit(1);
    }
	*/
    if ((req1 = (HECMW_Request *)calloc(pesize, sizeof(HECMW_Request)))
	== NULL) 
	   HECMW_vis_memory_exit("send_recv: req1");
    if ((req2 = (HECMW_Request *)calloc(pesize, sizeof(HECMW_Request)))
	== NULL) 
	HECMW_vis_memory_exit("send_recv: req2");
    nflag = 1;
  }
  /* SEND */
  for (neib = 0; neib < pesize; neib++) {
/*    if (neib != 0) istart = stack_export[neib - 1];
    else istart = 0;
*/
	  if(neib!=my_rank) {
         inum   = stack_export[neib+1] - stack_export[neib];
    
          for (k = stack_export[neib]; k < stack_export[neib] + inum; k++) {
              ws[k] = x[nod_export[k]];
		  }
          HECMW_Isend(&ws[stack_export[neib]], inum, HECMW_INT, neib, 0, repart_comm, &req1[neib]);
	  }
  }

  /* RECEIVE */
  for (neib = 0; neib < pesize; neib++) {
/*    if (neib != 0) istart = stack_import[neib-1];
    else istart = 0;
 */
   inum = stack_import[neib+1] - stack_import[neib];
   if(neib!=my_rank)
    HECMW_Irecv(&wr[stack_import[neib]], inum, HECMW_INT,
	      neib, 0, repart_comm, &req2[neib]);
  }
  HECMW_Barrier(repart_comm);
  /*
  HECMW_Waitall(neibpetot, req2, sta2);
  */
  for (neib = 0; neib < pesize; neib++) {
/*    if (neib != 0) istart = stack_import[neib-1];
    else istart = 0;
*/
    inum = stack_import[neib+1] - stack_import[neib];
	if(neib!=my_rank) {
    for (k = stack_import[neib]; k < stack_import[neib]+inum; k++) 
      y[nod_import[k]] = wr[k];
    }
	else {
		for(k=0;k<inum;k++)
/*		for (k = stack_import[neib]; k < stack_import[neib]+inum; k++)
*/
			y[nod_import[stack_import[neib]+k]]=x[nod_export[stack_export[neib]+k]];
	}
  }
  HECMW_Barrier(repart_comm);
  /*
  HECMW_Waitall(neibpetot, req1, sta1);
  
  free(sta1);
  free(sta2);
*/  
  free(req1);
  free(req2);
  free(ws);
  free(wr);
  return;
}


void int2_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, 
		     int *stack_export, 
		     int *x, int *y, 
		     HECMW_Comm repart_comm, int my_rank)
{
  /* Important:: node ID in nod_import and nod_export are all starting from 0  */

  HECMW_Status	*sta1, *sta2;
  HECMW_Request	*req1, *req2;

  int	nflag = 0;
  int	neib;
  int	istart, inum;
  int	k;
  int   *ws, *wr;

  if (nflag == 0) {
    sta1 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE*pesize, sizeof(HECMW_Status));
    if (sta1 == NULL) 
		HECMW_vis_memory_exit("send_recv: sta1");
    sta2 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE*pesize, sizeof(HECMW_Status));
    if (sta2 == NULL) 
		HECMW_vis_memory_exit("send_recv: sta2");
	
    if ((req1 = (HECMW_Request *)calloc(pesize, sizeof(HECMW_Request)))
	== NULL) 
	HECMW_vis_memory_exit("send_recv: req1");
    if ((req2 = (HECMW_Request *)calloc(pesize, sizeof(HECMW_Request)))
	== NULL) 
	HECMW_vis_memory_exit("send_recv: req2");
    nflag = 1;
  }
  /* SEND */
  for (neib = 0; neib < pesize; neib++) {
	  if(neib!=my_rank) {
         inum   = stack_export[neib+1] - stack_export[neib];
		 if(inum!=0) {
          HECMW_Isend(&x[stack_export[neib]], inum, HECMW_INT, neib, 0, repart_comm, &req1[neib]);
/*		  fprintf(stderr, "sending from PE %d to PE %d nodes %d\n", my_rank, neib, inum);
*/
		 }
	  }
  }

  /* RECEIVE */
  for (neib = 0; neib < pesize; neib++) {
/*    if (neib != 0) istart = stack_import[neib-1];
    else istart = 0;
 */
   inum = stack_import[neib+1] - stack_import[neib];
   if((neib!=my_rank) && (inum!=0)) {
    HECMW_Irecv(&y[stack_import[neib]], inum, HECMW_INT,
	      neib, 0, repart_comm, &req2[neib]);
/*	fprintf(stderr, "recv: PE %d from PE %d nodes %d\n", my_rank, neib, inum);
*/
   }
  }

  HECMW_Waitall(pesize, req2, sta2);


    inum = stack_import[my_rank+1] - stack_import[my_rank];
		for(k=0;k<inum;k++)
			y[stack_import[my_rank]+k]=x[stack_export[my_rank]+k];
 
  HECMW_Waitall(pesize, req1, sta1);
  
  free(sta1);
  free(sta2);
  
  free(req1);
  free(req2);
  return;
}

void double2_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, 
		     int *stack_export, double *x, double *y, 
		     HECMW_Comm repart_comm, int my_rank)
{
  /* Important:: node ID in nod_import and nod_export are all starting from 0  */

  HECMW_Status	*sta1, *sta2;
  HECMW_Request	*req1, *req2;

  int	nflag = 0;
  int	neib;
  int	istart, inum;
  int	k;

  if (nflag == 0) {
    sta1 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE*pesize, sizeof(HECMW_Status));
    if (sta1 == NULL) 
		HECMW_vis_memory_exit("send_recv: sta1");
    sta2 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE*pesize, sizeof(HECMW_Status));
    if (sta2 == NULL) 
		HECMW_vis_memory_exit("send_recv: sta2");
	
    if ((req1 = (HECMW_Request *)calloc(pesize, sizeof(HECMW_Request)))
	== NULL) 
	HECMW_vis_memory_exit("send_recv: req1");
    if ((req2 = (HECMW_Request *)calloc(pesize, sizeof(HECMW_Request)))
	== NULL) 
	HECMW_vis_memory_exit("send_recv: req2");
    nflag = 1;
  }
  /* SEND */
  for (neib = 0; neib < pesize; neib++) {
	  if(neib!=my_rank) {
         inum   = stack_export[neib+1] - stack_export[neib];
		 if(inum!=0) {
/*		  fprintf(stderr, "sending from PE %d to PE %d nodes %d\n", my_rank, neib, inum);
*/
          HECMW_Isend(&x[stack_export[neib]], inum, HECMW_DOUBLE, neib, 0, repart_comm, &req1[neib]);

		 }
	  }
  }

  /* RECEIVE */
  for (neib = 0; neib < pesize; neib++) {
/*    if (neib != 0) istart = stack_import[neib-1];
    else istart = 0;
 */
   inum = stack_import[neib+1] - stack_import[neib];
   if((neib!=my_rank) && (inum!=0)) {
/*	fprintf(stderr, "recv: PE %d from PE %d nodes %d\n", my_rank, neib, inum);
*/
    HECMW_Irecv(&y[stack_import[neib]], inum, HECMW_DOUBLE,
	      neib, 0, repart_comm, &req2[neib]);

   }
  }

  HECMW_Waitall(pesize, req2, sta2);


    inum = stack_import[my_rank+1] - stack_import[my_rank];
		for(k=0;k<inum;k++)
			y[stack_import[my_rank]+k]=x[stack_export[my_rank]+k];
 
  HECMW_Waitall(pesize, req1, sta1);
  
  free(sta1);
  free(sta2);
  
  free(req1);
  free(req2);
  return;
}

void int3_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, 
		     int *stack_export, 
		     int *x, int *y, 
		     HECMW_Comm repart_comm, int my_rank)
{
  /* Important:: node ID in nod_import and nod_export are all starting from 0  */

  HECMW_Status	*sta1, *sta2;
  HECMW_Request	*req1, *req2;
  HECMW_Status	stat;

  int	nflag = 0;
  int	neib;
  int	istart, inum;
  int	i, j,k;
  int   *ws, *wr;

	for(i=0;i<pesize;i++) {
		if(i!=my_rank) {
			inum=stack_export[i+1]-stack_export[i];
			if(inum!=0)
		   HECMW_Send(&x[stack_export[neib]], inum, HECMW_INT, i, 0, repart_comm);
		}
		else if(i==my_rank) {
         inum = stack_import[my_rank+1] - stack_import[my_rank];
		for(k=0;k<inum;k++)
			y[stack_import[my_rank]+k]=x[stack_export[my_rank]+k];
			for(j=0;j<pesize;j++) {
				if(j!=my_rank) {
					   inum = stack_import[j+1] - stack_import[j];
					   if(inum!=0)

					  HECMW_Recv(&y[stack_import[j]], inum, HECMW_INT, j, HECMW_ANY_TAG, repart_comm, &stat);
				}
			}
		}
	}

   HECMW_Barrier(repart_comm);
 
 
  return;
}

void double_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, int *nod_import,
		     int *stack_export, int *nod_export, 
		     double *x, double *y, 
		     HECMW_Comm repart_comm, int my_rank)
{
/*  HECMW_Status	*sta1, *sta2;
*/
  HECMW_Request	*req1, *req2;

  int	nflag = 0;
  int	neib;
  int	istart, inum;
  int	k;
  double   *ws, *wr;

  ws=(double *)calloc(n1, sizeof(int));
  wr=(double *)calloc(n2, sizeof(int));
  if((ws==NULL) || (wr==NULL)) 
	  HECMW_vis_memory_exit("send_recv: ws,wr");
  if (nflag == 0) {
/*    sta1 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status));
    if (sta1 == NULL) {
      fprintf(stderr, "Not enough memory\n");
      exit(1);
    }
    sta2 = (HECMW_Status *)calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status));
    if (sta2 == NULL) {
      fprintf(stderr, "Not enough memory\n");
      exit(1);
    }
	*/
    if ((req1 = (HECMW_Request *)calloc(pesize, sizeof(HECMW_Request)))
	== NULL) 
	HECMW_vis_memory_exit("send_recv: req1");
    if ((req2 = (HECMW_Request *)calloc(pesize, sizeof(HECMW_Request)))
	== NULL) 
	HECMW_vis_memory_exit("send_recv: req2");
    nflag = 1;
  }
  /* SEND */
  for (neib = 0; neib < pesize; neib++) {
/*    if (neib != 0) istart = stack_export[neib - 1];
    else istart = 0;
*/
	  if(neib!=my_rank) {
         inum   = stack_export[neib+1] - stack_export[neib];
    
          for (k = stack_export[neib]; k < stack_export[neib] + inum; k++) {
              ws[k] = x[nod_export[k]-1];
		  }
          HECMW_Isend(&ws[stack_export[neib]], inum, HECMW_DOUBLE, neib, 0, repart_comm, &req1[neib]);
	  }
  }

  /* RECEIVE */
  for (neib = 0; neib < pesize; neib++) {
/*    if (neib != 0) istart = stack_import[neib-1];
    else istart = 0;
 */

    inum = stack_import[neib+1] - stack_import[neib];
	if(neib!=my_rank) {
    for (k = stack_import[neib]; k < stack_import[neib]+inum; k++) 
      y[k] = wr[k];
    }
	else {
		for(k=0;k<inum;k++)
/*		for (k = stack_import[neib]; k < stack_import[neib]+inum; k++)
*/
			y[stack_import[neib]+k]=x[stack_export[neib]+k];
	}
  }

  HECMW_Barrier(repart_comm);
  /*
  HECMW_Waitall(neibpetot, req1, sta1);
  
  free(sta1);
  free(sta2);
*/  
  free(req1);
  free(req2);
  free(ws);
  free(wr);
  return;
}
