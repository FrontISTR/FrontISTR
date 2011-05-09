//
//  ContactNode.cpp
//
//
//
//              2009.10.26
//              2009.10.26
//              k.Takeda
#include "ContactNode.h"
using namespace pmw;


// construct & destruct
//
CContactNode::CContactNode()
{
    mbMesh=false;
    mbNode=false;
}
CContactNode::~CContactNode()
{
    ;
}

//// 節点集合の"SkinFaceID",とFace内部の"ノード配列インデックス番号"をセット
////
//void CContactNode::pushSkinFaceID(const uint& skinFaceID, const uint& nodeIndex)
//{
//    mvAggFaceID.push_back(skinFaceID);
//    mmFaceVertNum[skinFaceID]= nodeIndex;
//}


// EQUATION 関連パラメータ
// ------------
// Arbitrary DOF(任意自由度)
//
// "変位"の配列確保
//
void CContactNode::resizeDisp(const uint& dof)
{
    mvDisplacement.resize(dof);
}

// 変位の初期化
//
void CContactNode::initDisp()
{
    uint numOfDOF;
    uint idof;
    
    numOfDOF= mvDisplacement.size();
    for(idof=0; idof< numOfDOF; idof++){
        mvDisplacement[idof]=0.0;
    };
}

// スカラーパラメータの配列確保(複数のスカラー)
//
void CContactNode::resizeScalar(const uint& numOfScalar)
{
    mvScalar.resize(numOfScalar);
}

// スカラーの初期化
//
void CContactNode::initScalar()
{
    uint numOfScalar= mvScalar.size();
    uint i;

    for(i=0; i< numOfScalar; i++){
        mvScalar[i]= 0.0;
    };
}




