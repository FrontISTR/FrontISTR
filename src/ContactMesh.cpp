//
// ContactMesh.cpp
//
//			2009.01.08
//			2009.01.08
//			k.Takeda
#include "ContactMesh.h"
using namespace pmw;

//
//
CContactMesh::CContactMesh()
{
}

CContactMesh::~CContactMesh()
{
  for_each(mvContactElement.begin(), mvContactElement.end(), DeleteObject());

  //debug
  cout << "~CContactMesh" << endl;
}
