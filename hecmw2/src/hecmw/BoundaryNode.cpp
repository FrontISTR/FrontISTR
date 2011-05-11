//
//  BoundaryNode.cpp
//
//  境界Node
//
//                          2009.05.13
//                          2009.05.13
//                          k.Takeda
#include <vector>

#include "BoundaryNode.h"
using namespace pmw;

CBoundaryNode::CBoundaryNode()
{
    ;
}
CBoundaryNode::~CBoundaryNode()
{
//    //debug
//    cout << "~CBoundaryNode" << endl;
}


// Level数分の配列確保
// ----
void CBoundaryNode::resizeValue(const uiint& numOfDiffLevel)
{
    mvValue.resize(numOfDiffLevel);
}


// ゼロクリア 
// ----
void CBoundaryNode::initValue(const uiint& dof, const uiint& mgLevel)
{
    uiint diffLevel= mgLevel - mMGLevel;//diffLevelは,計算Levelと自身の差
    
    mvValue[diffLevel][dof]= 0.0;
}


// 代入
// ---
void CBoundaryNode::setValue(const uiint& dof, const uiint& mgLevel, const double& val)
{
    uiint diffLevel= mgLevel - mMGLevel;

    //cout << "BoundaryNode::setValue, diffLevel= " << diffLevel << ", mvValue.size=" << mvValue.size() << endl;

    mvValue[diffLevel][dof]= val;
}


// 加算(足し合わせ)
// ----
void CBoundaryNode::addValue(const uiint& dof, const uiint& mgLevel, const double& val)
{
    uiint diffLevel= mgLevel - mMGLevel;

    mvValue[diffLevel][dof] += val;
}


// 境界値の提供
// ----
double& CBoundaryNode::getValue(const uiint& dof, const uiint& mgLevel)
{
    uiint diffLevel= mgLevel - mMGLevel;

    return mvValue[diffLevel][dof];
}



