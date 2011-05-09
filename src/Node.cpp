//
// Node3.cpp
//
//				2008.12.18
//				2008.12.18
//				k.Takeda
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


