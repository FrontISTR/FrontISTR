/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Material.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
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
