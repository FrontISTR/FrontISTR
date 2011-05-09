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
void CBoundaryNode::resizeValue(const uint& numOfDiffLevel)
{
    mvValue.resize(numOfDiffLevel);
}


// ゼロクリア 
// ----
void CBoundaryNode::initValue(const uint& dof, const uint& mgLevel)
{
    uint diffLevel= mgLevel - mMGLevel;//diffLevelは,計算Levelと自身の差
    
    mvValue[diffLevel][dof]= 0.0;
}


// 代入
// ---
void CBoundaryNode::setValue(const uint& dof, const uint& mgLevel, const double& val)
{
    uint diffLevel= mgLevel - mMGLevel;

    //cout << "BoundaryNode::setValue, diffLevel= " << diffLevel << ", mvValue.size=" << mvValue.size() << endl;

    mvValue[diffLevel][dof]= val;
}


// 加算(足し合わせ)
// ----
void CBoundaryNode::addValue(const uint& dof, const uint& mgLevel, const double& val)
{
    uint diffLevel= mgLevel - mMGLevel;

    mvValue[diffLevel][dof] += val;
}


// 境界値の提供
// ----
double& CBoundaryNode::getValue(const uint& dof, const uint& mgLevel)
{
    uint diffLevel= mgLevel - mMGLevel;

    return mvValue[diffLevel][dof];
}



