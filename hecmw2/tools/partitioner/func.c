/*
 ----------------------------------------------------------
|
| Software Name :part Ver 0.1 beta
|
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include <stdlib.h>
#include "func.h"

int getCommMeshID(int myrank, int rank, int n)
{
    //printf("getCommMeshID  ");
    int x,y;
    if(myrank < rank) {
        x = rank;
        y = myrank;
    } else if(rank < myrank) {
        x = myrank;
        y = rank;
    } else {
        exit(-1);
    }

    int id = 0;
    int i,j;
    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
            if(i == y && j == x) {
                break;
            }
            if(j > i) {
                id++;
            }
        }
        if(i == y && j == x) {
            break;
        }
    }
    //printf("%d %d %d\n", myrank, rank, id);
    return id;
}
