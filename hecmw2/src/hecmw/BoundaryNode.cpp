/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BoundaryNode.cpp
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
#include <vector>
#include "BoundaryNode.h"
using namespace pmw;
CBoundaryNode::CBoundaryNode()
{
    ;
}
CBoundaryNode::~CBoundaryNode()
{
}
void CBoundaryNode::resizeValue(const uiint& numOfDiffLevel)
{
    mvValue.resize(numOfDiffLevel);
}
void CBoundaryNode::initValue(const uiint& dof, const uiint& mgLevel)
{
    uiint diffLevel= mgLevel - mMGLevel;
    mvValue[diffLevel][dof]= 0.0;
}
void CBoundaryNode::setValue(const uiint& dof, const uiint& mgLevel, const double& val)
{
    uiint diffLevel= mgLevel - mMGLevel;
    mvValue[diffLevel][dof]= val;
}
void CBoundaryNode::addValue(const uiint& dof, const uiint& mgLevel, const double& val)
{
    uiint diffLevel= mgLevel - mMGLevel;
    mvValue[diffLevel][dof] += val;
}
double& CBoundaryNode::getValue(const uiint& dof, const uiint& mgLevel)
{
    uiint diffLevel= mgLevel - mMGLevel;
    return mvValue[diffLevel][dof];
}
