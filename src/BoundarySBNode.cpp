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
void CBoundarySBNode::addDOF(const uint& dof)
{
    mvDOF.push_back(dof);
}
uint& CBoundarySBNode::getDOF(const uint& index)
{
    return mvDOF[index];
}
uint CBoundarySBNode::getNumOfDOF()
{
    return mvDOF.size();
}

// 境界値
// ----
void CBoundarySBNode::setValue(const uint& dof, const double& val)
{
    mmValue[dof]= val;
}
double& CBoundarySBNode::getValue(const uint& dof)
{
    return mmValue[dof];
}






