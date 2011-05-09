//
//  AssyModel.cpp
//
//			2009.04.20
//			2008.11.05
//			k.Takeda
#include "AssyModel.h"
using namespace pmw;

CAssyModel::CAssyModel(void)
{
//    mvContactMesh.resize(0);
//    mvMesh.resize(0);
}

CAssyModel::~CAssyModel(void)
{
    // Mesh, ContactMesh
    for_each(mvMesh.begin(), mvMesh.end(), DeleteObject());
    for_each(mvContactMesh.begin(), mvContactMesh.end(), DeleteObject());

    //debug
    cout << "~CAssyModel" << endl;
}



