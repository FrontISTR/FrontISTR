/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Matrix.cpp
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
