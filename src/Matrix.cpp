/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Matrix.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
/*
 * Matrix.cpp
 *
 *  Created on: Jul 23, 2009
 *      Author: goto
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
