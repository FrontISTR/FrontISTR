/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Hexa2.cpp
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
#include "Hexa2.h"
using namespace pmw;
uiint CHexa2::mnElemType = ElementType::Hexa2;
uiint CHexa2::mnElemOrder = 2;
uiint CHexa2::mNumOfFace =  6;
uiint CHexa2::mNumOfEdge = 12;
uiint CHexa2::mNumOfNode = 20;
uiint CHexa2::mNumOfVert =  8;
CHexa2::CHexa2()
{
    ;
}
CHexa2::~CHexa2()
{
    ;
}
void CHexa2::initialize()
{
    mvNode.resize(mNumOfNode);
    mvvEdgeElement.resize(mNumOfEdge);
    mvEdgeInterNode.resize(mNumOfEdge);
    mvb_edge = new bool[mNumOfEdge];
    uiint i;
    for(i=0; i< mNumOfEdge; i++) {
        mvb_edge[i]=false;
    };
    mvFaceElement.resize(mNumOfFace);
    mvFaceNode.resize(mNumOfFace);
    mvb_face = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvb_face[i]=false;
    };
    mvProgElement.resize(mNumOfNode);
    mvbMPCFace = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvbMPCFace[i]= false;
    };
////    mvbCommEntity = new bool[mNumOfFace];
////    for(i=0; i< mNumOfFace; i++){
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
void CHexa2::replaseEdgeNode()
{
    CNode *pNode;
    uiint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++) {
        pNode = mvEdgeInterNode[iedge];
        mvNode[mNumOfVert + iedge] = pNode;
        mmIDLocal[pNode->getID()] = mNumOfVert + iedge;
    };
}
