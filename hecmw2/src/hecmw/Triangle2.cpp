/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Triangle2.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include <vector>
#include "Triangle2.h"
using namespace pmw;
uiint CTriangle2::mnElemType = ElementType::Triangle2;
uiint CTriangle2::mnElemOrder = 2;
uiint CTriangle2::mNumOfFace = 1;
uiint CTriangle2::mNumOfEdge = 3;
uiint CTriangle2::mNumOfNode = 6;
uiint CTriangle2::mNumOfVert = 3;
CTriangle2::CTriangle2()
{
    ;
}
CTriangle2::~CTriangle2()
{
    ;
}
void CTriangle2::initialize()
{
    mvNode.resize(mNumOfNode);
    mvvEdgeElement.resize(mNumOfEdge);
    mvEdgeInterNode.resize(mNumOfEdge);
    mvb_edge = new bool[mNumOfEdge];
    uiint i;
    for(i=0; i< mNumOfEdge; i++) {
        mvb_edge[i] = false;
    };
    mvFaceElement.resize(mNumOfFace);
    mvFaceNode.resize(mNumOfFace);
    mvb_face = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvb_face[i] = false;
    };
    mvProgElement.resize(mNumOfNode);
    mvbMPCFace = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvbMPCFace[i]= false;
    };
////    mvbCommEntity = new bool[mNumOfEdge];
////    for(i=0; i< mNumOfEdge; i++){
////        mvbCommEntity[i]= false;
////    };
    mvbCommLFace = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvbCommLFace[i]= false;
    };
    mvbCommLEdge = new bool[mNumOfEdge];
    for(i=0; i< mNumOfEdge; i++) {
        mvbCommLEdge[i]= false;
    };
    mvbCommLVert = new bool[mNumOfVert];
    for(i=0; i< mNumOfVert; i++) {
        mvbCommLVert[i]= false;
    }
}
void CTriangle2::replaseEdgeNode()
{
    uiint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++) {
        mvNode[mNumOfVert + iedge] = mvEdgeInterNode[iedge];
    };
}
