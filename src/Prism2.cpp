//
// Prism2.cpp
//
//              2010.11.19
//              k.Takeda
#include <vector>

#include "Prism2.h"
using namespace pmw;


uint CPrism2::mnElemType =  ElementType::Prism2;
uint CPrism2::mnElemOrder = 2;
uint CPrism2::mNumOfFace =  5;
uint CPrism2::mNumOfEdge =  9;
uint CPrism2::mNumOfNode = 15;
uint CPrism2::mNumOfVert =  6;

CPrism2::CPrism2()
{
    ;
}
CPrism2::~CPrism2()
{
    ;
}
void CPrism2::initialize()
{
    mvNode.resize(mNumOfNode);

    mvvEdgeElement.resize(mNumOfEdge);

    mvEdgeInterNode.resize(mNumOfEdge);

    mvb_edge = new bool[mNumOfEdge];
    uint i;
    for(i=0; i< mNumOfEdge; i++){
      mvb_edge[i] = false;
    };

    mvFaceElement.resize(mNumOfFace);
    mvFaceNode.resize(mNumOfFace);

    mvb_face = new bool[mNumOfFace];
    // Face ノード bool
    //
    for(i=0; i< mNumOfFace; i++){
        mvb_face[i] = false;
    };

    //CommElementのprolongation用
    mvProgElement.resize(mNumOfNode);

    //MPC属性
    mvbMPCFace = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++){
        mvbMPCFace[i]= false;
    };

    //通信界面(CommMesh2)属性
    mvbCommEntity = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++){
        mvbCommEntity[i]= false;
    };
}

//
// 辺ノードを、2次ノードへ入れ替え
//
void CPrism2::replaseEdgeNode()
{
    CNode *pNode;
    uint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        pNode = mvEdgeInterNode[iedge];
        mvNode[mNumOfVert + iedge] = pNode;

        mmIDLocal[pNode->getID()] = mNumOfVert + iedge;
    };

    //ReInit_IDLocal();
}

////
//// 1. 辺-面 Element*配列を解放
//// 2. 辺-面 Node*   配列を解放 (2次要素は辺ノードを残す)
////
//void CPrism2::deleteProgData()
//{
//    // Edge
//    uint iedge;
//    for(iedge=0; iedge < mNumOfEdge; iedge++){
//        vector<CElement*>().swap(mvvEdgeElement[iedge]);
//    };
//    vector<vector<CElement*> >().swap(mvvEdgeElement);
//    //vector<CNode*>().swap(mvEdgeInterNode);//2次要素の場合は残すこと
//    delete []mvb_edge;
//
//
//    // Face
//    vector<CElement*>().swap(mvFaceElement);
//    vector<CNode*>().swap(mvFaceNode);
//    delete []mvb_face;
//
//
//    delete []mvbMPCFace;   // 面番号のどれがMPCするのか
//    delete []mvbCommEntity;// CommMesh2の属性
//}





