/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Node.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Node.h"

CNode::CNode(void)
{
	for(int i=0; i < 3; i++) mdCoord[i]= 0.0;
}

CNode::~CNode(void)
{
	;
}
void CNode::setID(size_t id)
{
	mID= id;
}
void CNode::setCoord(const double& x, const double& y, const double& z)
{
	mdCoord[0]=x;  mdCoord[1]=y;  mdCoord[2]=z;
}


