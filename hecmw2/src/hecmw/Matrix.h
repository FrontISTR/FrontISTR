/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Matrix.h
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
