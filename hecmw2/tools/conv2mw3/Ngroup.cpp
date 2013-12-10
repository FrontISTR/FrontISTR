/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Ngroup.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Ngroup.h"

CNgroup::CNgroup()
{
}
CNgroup::~CNgroup()
{
}
//--
// 
//--
void CNgroup::addNode(CNode* pNode)
{
	mvNode.push_back(pNode);
}
size_t CNgroup::getNumOfNode()
{
	return mvNode.size();
}
CNode* CNgroup::getNode(size_t index)
{
	return mvNode[index];
}

