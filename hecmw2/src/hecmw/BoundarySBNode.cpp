//
//  BoundarySBNode.cpp
//
//
//          2010.06.25
//          k.Takeda
#include "BoundarySBNode.h"
using namespace pmw;


CBoundarySBNode::CBoundarySBNode()
{
    ;
}
CBoundarySBNode::~CBoundarySBNode()
{
    ;
}

// DOF
// ---
void CBoundarySBNode::addDOF(const uiint& dof)
{
    mvDOF.push_back(dof);
}
uiint& CBoundarySBNode::getDOF(const uiint& index)
{
    return mvDOF[index];
}
uiint CBoundarySBNode::getNumOfDOF()
{
    return mvDOF.size();
}

// 境界値
// ----
void CBoundarySBNode::setValue(const uiint& dof, const double& val)
{
    mmValue[dof]= val;
}
double& CBoundarySBNode::getValue(const uiint& dof)
{
    return mmValue[dof];
}






