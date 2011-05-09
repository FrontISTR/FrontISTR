//
// DiscreteVolume.cpp
//
//
//          2010.06.22
//          k.Takeda
#include "DiscreteVolume.h"
using namespace pmw;

CDiscreteVolume::CDiscreteVolume()
{
    // Prismの頂点番号で表した3個のTetra
    uint prismDiscre[3][4]={
        {0,4,3,5},
        {0,2,4,5},
        {0,2,1,4}
    };
    // Hexaの頂点番号で表した6個のTetra
    uint hexaDiscre[6][4]={
        {0,7,4,5},
        {0,7,5,1},
        {0,3,7,1},
        {7,3,6,5},
        {3,6,5,1},
        {1,2,3,6}
    };

    // classメンバーへコピー
    uint i,ii;
    for(i=0; i< 3; i++){
        for(ii=0; ii< 4; ii++){
            mPrismDiscre[i][ii]= prismDiscre[i][ii];
        };
    };
    for(i=0; i< 6; i++){
        for(ii=0; ii< 4; ii++){
            mHexaDiscre[i][ii]= hexaDiscre[i][ii];
        };
    };
}
CDiscreteVolume::~CDiscreteVolume()
{
    ;
}










