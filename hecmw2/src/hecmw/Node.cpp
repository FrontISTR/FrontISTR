/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Node.cpp
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
#include "Node.h"
using namespace pmw;
CNode::CNode(void)
{
    mbSComm = false;
}
CNode::~CNode(void)
{
    ;
}
void CNode::markingSCommNode()
{
    mbSComm = true;
}
bool CNode::isSCommNode()
{
    return mbSComm;
}
