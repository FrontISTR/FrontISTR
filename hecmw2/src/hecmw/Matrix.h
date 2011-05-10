/*
 * Matrix.h
 *
 *  Created on: Jul 23, 2009
 *      Author: goto
 */

#ifndef MATRIX_H_
#define MATRIX_H_

namespace pmw
{

class CVector;

class CMatrix
{
public:
    CMatrix();
    virtual ~CMatrix();

    virtual void multVector(CVector *pV, CVector *pP) const = 0;
    void residual(CVector *pV, const CVector *pF, CVector *pR) const;

private:
};

}

#endif /* MATRIX_H_ */
