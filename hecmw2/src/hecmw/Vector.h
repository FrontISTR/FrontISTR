/*
 * Vector.h
 *
 *  Created on: Jul 23, 2009
 *      Author: goto
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
        //////////////////////////////////////
	typedef ublas::vector<double> ElemType;
        //////////////////////////////////////

	CVector(CMesh *pMesh, const uiint& nDOF);
	CVector(const CVector *pVector);
	virtual ~CVector();

	uiint size() const;
	// size_t lenInternal() const;
	const ElemType &operator[](uiint idx) const;
	ElemType &operator[](uiint idx);

        void Vector_Clear();// Matrix 0 clear

	void setZero();
	void setValue(uiint inode, uiint idof, double value);
	void addValue(uiint inode, uiint idof, double value);
	double& getValue(uiint inode, uiint idof);
	void sumSV(double alpha, const CVector *pX, CVector *pY) const;
	void addSV(double alpha, const CVector *pX);
	void add(const CVector *pX);
	void subst(const CVector *pX);
//	CMesh getMesh(){return mpMesh;};

	double norm2() const;
	double innerProd(const CVector *pX) const;
	void updateCommBoundary();

	uiint restrictTo(CVector *pV) const;
	uiint prolongateFrom(const CVector *pV);

	void print_elem() const ;//{
	//	for(uint i=0;i<size();i++) printf("%d %e %e %e \n",i,mvVector[i](0),mvVector[i](1),mvVector[i](2));
	//};

        void dump();//2011.01.12 列ベクトルのダンプ

private:
	uiint mnNode;
	uiint mnNodeInternal;
	int mnDOF;
	std::vector<ElemType> mvVector;
	// restrictor
	// prolongator
	CMesh *mpMesh;

        bool isScopeNode(const uiint& idx) const;//Vector範囲内に存在するか.
};

}

#endif /* VECTOR_H_ */
