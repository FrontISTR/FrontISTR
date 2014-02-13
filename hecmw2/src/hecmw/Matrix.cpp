/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Matrix.cpp
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
#include "Matrix.h"
#include "Vector.h"
namespace pmw
{
CMatrix::CMatrix()
{
}
CMatrix::~CMatrix()
{
}
void CMatrix::residual(CVector *pV, const CVector *pF, CVector *pR) const
{
    multVector(pV, pR);
    int len = pV->size();
    for (int i = 0; i < len; i++) {
        (*pR)[i] = (*pF)[i] - (*pR)[i];
    }
}
}
