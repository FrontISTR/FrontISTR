/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BoundarySBNode.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
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
void CBoundarySBNode::setValue(const uiint& dof, const double& val)
{
    mmValue[dof]= val;
}
double& CBoundarySBNode::getValue(const uiint& dof)
{
    return mmValue[dof];
}
