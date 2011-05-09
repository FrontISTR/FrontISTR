/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Material.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
void CMaterial::setValue(const uint& prop_Type, const double& value)
{
    mmValue[prop_Type]= value;
}
