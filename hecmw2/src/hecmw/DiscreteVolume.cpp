/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/DiscreteVolume.cpp
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
#include "DiscreteVolume.h"
using namespace pmw;
CDiscreteVolume::CDiscreteVolume()
{
    uiint prismDiscre[3][4]= {
        {0,4,3,5},
        {0,2,4,5},
        {0,2,1,4}
    };
    uiint hexaDiscre[6][4]= {
        {0,7,4,5},
        {0,7,5,1},
        {0,3,7,1},
        {7,3,6,5},
        {3,6,5,1},
        {1,2,3,6}
    };
    uiint i,ii;
    for(i=0; i< 3; i++) {
        for(ii=0; ii< 4; ii++) {
            mPrismDiscre[i][ii]= prismDiscre[i][ii];
        };
    };
    for(i=0; i< 6; i++) {
        for(ii=0; ii< 4; ii++) {
            mHexaDiscre[i][ii]= hexaDiscre[i][ii];
        };
    };
}
CDiscreteVolume::~CDiscreteVolume()
{
    ;
}
