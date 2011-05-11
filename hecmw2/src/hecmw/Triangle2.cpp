//
// Triangle2.cpp
//
//              2010.11.19
//              k.Takeda
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

    //通信界面(CommMesh2)属性:辺
    mvbCommEntity = new bool[mNumOfEdge];
    for(i=0; i< mNumOfEdge; i++){
        mvbCommEntity[i]= false;
    };
}

//
// 辺ノードを、2次ノードへ入れ替え
//
void CTriangle2::replaseEdgeNode()
{
    uiint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        mvNode[mNumOfVert + iedge] = mvEdgeInterNode[iedge];
    };
}

////
//// 1. 辺-面 Element*配列を解放
//// 2. 辺-面 Node*   配列を解放 (2次要素は辺ノードを残す)
////
//void CTriangle2::deleteProgData()
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







