//

#include <vector>

// Beam2.cpp
//
//              2010.11.19
//              k.Takeda
#include "Beam2.h"
#include "AggregateElement.h"
using namespace pmw;

uint CBeam2::mnElemType = ElementType::Beam2;
uint CBeam2::mnElemOrder = 2;
uint CBeam2::mNumOfFace = 0;
uint CBeam2::mNumOfEdge = 1;
uint CBeam2::mNumOfNode = 3;
uint CBeam2::mNumOfVert = 2;

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
    uint i;
    for(i=0; i< mNumOfEdge; i++){
      mvb_edge[i] = false;
    };

    //CommElementのprolongation用
    mvProgElement.resize(mNumOfNode);


    //通信界面(CommMesh2)属性:点
    mvbCommEntity = new bool[mNumOfNode];
    for(i=0; i< mNumOfNode; i++){
        mvbCommEntity[i]= false;
    };
}

//
// 辺ノードを、2次ノードへ入れ替え
//
void CBeam2::replaseEdgeNode()
{
    uint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        mvNode[mNumOfVert + iedge] = mvEdgeInterNode[iedge];
    };
}

////
//// 1. 辺-面 Element*配列を解放
//// 2. 辺-面 Node*   配列を解放 (2次要素は辺ノードを残す)
////
//void CBeam2::deleteProgData()
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
//    delete []mvbCommEntity;
//}















