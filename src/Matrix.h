/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Matrix.h
|
|                     Written by T.Takeda,    2010/06/01
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
