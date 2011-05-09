/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Vector.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef VECTOR_H_
#define VECTOR_H_
#include "TypeDef.h"
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
using namespace boost::numeric;
namespace pmw
{
class CMesh;
class CVector
{
public:
	typedef ublas::vector<double> ElemType;
	CVector(/* const */ CMesh *pMesh);
	CVector(const CVector *pVector);
	virtual ~CVector();
	size_t size() const;
	const ElemType &operator[](size_t idx) const;
	ElemType &operator[](size_t idx);
	void setZero();
	void setValue(int inode, int idof, double value);
	void addValue(int inode, int idof, double value);
	double getValue(int inode, int idof);
	void sumSV(double alpha, const CVector *pX, CVector *pY) const;
	void addSV(double alpha, const CVector *pX);
	void add(const CVector *pX);
	void subst(const CVector *pX);
	double norm2() const;
	double innerProd(const CVector *pX) const;
	void updateCommBoundary();
	int restrictTo(CVector *pV) const;
	int prolongateFrom(const CVector *pV);
	void print_elem() const ;
private:
	size_t mnNode;
	size_t mnNodeInternal;
	int mnDOF;
	std::vector<ElemType> mvVector;
	CMesh *mpMesh;
};
}
#endif /* VECTOR_H_ */
