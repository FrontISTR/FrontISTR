/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Material.cpp
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
#include "Material.h"
using namespace pmw;
#include <iostream>
CMaterial::CMaterial()
{
    ;
}
CMaterial::~CMaterial()
{
    cout << "~CMaterial" << endl;
}
void CMaterial::setValue(const uiint& prop_Type, const double& value)
{
    mmValue[prop_Type]= value;
}
