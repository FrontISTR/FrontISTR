/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Beam2.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include <vector>
#include "Beam2.h"
#include "AggregateElement.h"
using namespace pmw;
uiint CBeam2::mnElemType = ElementType::Beam2;
uiint CBeam2::mnElemOrder = 2;
uiint CBeam2::mNumOfFace = 0;
uiint CBeam2::mNumOfEdge = 1;
uiint CBeam2::mNumOfNode = 3;
uiint CBeam2::mNumOfVert = 2;
CBeam2::CBeam2()
{
    ;
}
CBeam2::~CBeam2()
{
    ;
}
void CBeam2::initialize()
{
    mvNode.resize(mNumOfNode);
    mvvEdgeElement.resize(mNumOfEdge);
    mvEdgeInterNode.resize(mNumOfEdge);
    mvb_edge = new bool[mNumOfEdge];
    uiint i;
    for(i=0; i< mNumOfEdge; i++){
      mvb_edge[i] = false;
    };
    mvProgElement.resize(mNumOfNode);
    mvbCommEntity = new bool[mNumOfNode];
    for(i=0; i< mNumOfNode; i++){
        mvbCommEntity[i]= false;
    };
}
void CBeam2::replaseEdgeNode()
{
    uiint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        mvNode[mNumOfVert + iedge] = mvEdgeInterNode[iedge];
    };
}
