/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundarySBNode.cpp
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
//--
// 境界値
//--
void CBoundarySBNode::setValue(const uiint& dof, const double& val)
{
    mmValue[dof]= val;
}
double& CBoundarySBNode::getValue(const uiint& dof)
{
    return mmValue[dof];
}
//--
// 基礎データ(ディレクレ数式処理 基礎データ)
//--
void CBoundarySBNode::setEntValue(const uiint& dof, const double& val)
{
    mmEntValue[dof]= val;
}
double& CBoundarySBNode::getEntValue(const uiint& dof)
{
    return mmEntValue[dof];
}

double& CBoundarySBNode::getX()
{
    return mpNode->getX();
}
double& CBoundarySBNode::getY()
{
    return mpNode->getY();
}
double& CBoundarySBNode::getZ()
{
    return mpNode->getZ();
}

