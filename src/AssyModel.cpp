/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   AssyModel.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "AssyModel.h"
using namespace pmw;
CAssyModel::CAssyModel(void)
{
}
CAssyModel::~CAssyModel(void)
{
    for_each(mvMesh.begin(), mvMesh.end(), DeleteObject());
    for_each(mvContactMesh.begin(), mvContactMesh.end(), DeleteObject());
    cout << "~CAssyModel" << endl;
}
CMesh* CAssyModel::getMesh_ID(const uint& id)
{
    uint index= moBucketMesh.getIndexMesh(id);
    return mvMesh[index];
}
void CAssyModel::addContactMesh(CContactMesh* pContactMesh, const uint& id)
{
    mvContactMesh.push_back(pContactMesh);
    uint index= mvContactMesh.size()-1;
    mmContactID2Index[id]= index;
}
