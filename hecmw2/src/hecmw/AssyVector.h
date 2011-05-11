/*
 * AssyVector.h
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
 */

#ifndef ASSYVECTOR_H_
#define ASSYVECTOR_H_

#include <vector>
#include "Vector.h"

namespace pmw
{

class CAssyModel;

class CAssyVector
{
public:
	CAssyVector(CAssyModel *pAssyModel, const uiint& nDOF);
	CAssyVector(const CAssyVector *pAssyVector);
	virtual ~CAssyVector();

	uiint size() const;
	// int lenInternal() const;
	const CVector::ElemType &operator[](uiint idx) const;
	CVector::ElemType &operator[](uiint idx);

	const double &operator()(uiint meshID, uiint nodeID, uiint dof) const;
	double &operator()(uiint meshID, uiint nodeID, uiint dof);

        void Vector_Clear(const uiint& iMesh);// Matrix 0 clear

	void setZero();//use GMRES
	void setValue(uiint imesh, uiint inode, uiint idof, double value);
        void addValue(const uiint& imesh, const uiint& inode, const uiint& idof, const double& value);
	double& getValue(uiint imesh, uiint inode, uiint idof);
	void add(const CAssyVector *pV);
	void subst(const CAssyVector *pV);

	double norm2() const;
	double innerProd(const CAssyVector *pX) const;
	void updateCommBoundary();
	void sumupCommBoundary();

	uiint restrictTo(CAssyVector *pVc) const;
	uiint prolongateFrom(const CAssyVector *pVc);

	uiint getNumOfVector() const;
	const CVector *getVector(uiint index) const;
	CVector *getVector(uiint index);
	
	void print_elem() const {
		for (std::vector<CVector*>::const_iterator iv = mvVector.begin(); iv != mvVector.end(); iv++) (*iv)->print_elem();
	};

        uiint& getDOF();

        void dump();//2011.01.12 列ベクトルのダンプ

private:
	CAssyModel *mpAssyModel;

	std::vector<CVector*> mvVector;
	// commTable

        uiint mnDOF;//アセンブル方程式のDOF
};

}

#endif /* ASSYVECTOR_H_ */
