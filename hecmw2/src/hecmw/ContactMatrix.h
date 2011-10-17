/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/ContactMatrix.h
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
#ifndef CONTACTMATRIX_H_
#define CONTACTMATRIX_H_
namespace pmw
{
class CAssyVector;
class CContactMatrix
{
public:
	CContactMatrix();
	virtual ~CContactMatrix();
	void multVectorAdd(CAssyVector *pP, CAssyVector *pX);
};
}
#endif /* CONTACTMATRIX_H_ */
