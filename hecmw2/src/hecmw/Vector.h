/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Vector.h
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
	CVector(CMesh *pMesh, const uiint& nDOF);
	CVector(const CVector *pVector);
	virtual ~CVector();
	uiint size() const;
	const ElemType &operator[](uiint idx) const;
	ElemType &operator[](uiint idx);
        void Vector_Clear();
	void setZero();
	void setValue(uiint inode, uiint idof, double value);
	void addValue(uiint inode, uiint idof, double value);
	double& getValue(uiint inode, uiint idof);
	void sumSV(double alpha, const CVector *pX, CVector *pY) const;
	void addSV(double alpha, const CVector *pX);
	void add(const CVector *pX);
	void subst(const CVector *pX);
	double norm2() const;
	double innerProd(const CVector *pX) const;
	void updateCommBoundary();
	uiint restrictTo(CVector *pV) const;
	uiint prolongateFrom(const CVector *pV);
	void print_elem() const ;
        void dump();
private:
	uiint mnNode;
	uiint mnNodeInternal;
	int mnDOF;
	std::vector<ElemType> mvVector;
	CMesh *mpMesh;
        bool isScopeNode(const uiint& idx) const;
};
}
#endif /* VECTOR_H_ */
