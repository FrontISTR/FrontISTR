//
//  Material.cpp
//
//  材質データ
//
//                          2009.05.18
//                          2009.05.18
//                          k.Takeda
#include "Material.h"
using namespace pmw;

#include <iostream>
//
//
CMaterial::CMaterial()
{
    ;
}

CMaterial::~CMaterial()
{
    //debug
    cout << "~CMaterial" << endl;
}

// 材質データのセット
// --
void CMaterial::setValue(const uint& prop_Type, const double& value)
{
    mmValue[prop_Type]= value;
}



