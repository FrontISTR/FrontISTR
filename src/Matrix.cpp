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
    // TODO Auto-generated constructor stub

}

CMatrix::~CMatrix()
{
    // TODO Auto-generated destructor stub
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
