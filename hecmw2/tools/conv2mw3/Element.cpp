/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Element.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/

#include "Element.h"


CElement::CElement()
{
	;
}
CElement::~CElement()
{
	;
}
void CElement::setID(size_t id)
{
	mID= id;
}
void CElement::setNode(size_t index, CNode* pNode)
{
	mvNode[index]= pNode;
}



