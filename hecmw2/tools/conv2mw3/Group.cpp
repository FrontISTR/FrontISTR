/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Group.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Group.h"

CGroup::CGroup()
{
	mnBndType = BndType::NotUse;
}
CGroup::~CGroup()
{
}
void CGroup::setGroupName(string name)
{
	mGroupName= name;
}
void CGroup::setPartsName(string name)
{
	mPartsName= name;
}

void CGroup::setID(size_t id)
{
	mID=id;
}
size_t CGroup::getID()
{
	return mID;
}

void CGroup::setBndType(size_t nType)
{
	mnBndType= nType;
}
size_t CGroup::getBndType()
{
	return mnBndType;
}


