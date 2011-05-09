/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ContactMatrix.h
|
|                     Written by T.Takeda,    2010/06/01
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
