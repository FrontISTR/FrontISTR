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
	CAssyVector(CAssyModel *pAssyModel, const uint& nDOF);
	CAssyVector(const CAssyVector *pAssyVector);
	virtual ~CAssyVector();

	size_t size() const;
	// int lenInternal() const;
	const CVector::ElemType &operator[](size_t idx) const;
	CVector::ElemType &operator[](size_t idx);

	const double &operator()(size_t meshID, size_t nodeID, size_t dof) const;
	double &operator()(size_t meshID, size_t nodeID, size_t dof);

        void Vector_Clear(const uint& iMesh);// Matrix 0 clear

	void setZero();//use GMRES
	void setValue(int imesh, int inode, int idof, double value);
        void addValue(const uint& imesh, const uint& inode, const uint& idof, const double& value);
	double getValue(int imesh, int inode, int idof);
	void add(const CAssyVector *pV);
	void subst(const CAssyVector *pV);

	double norm2() const;
	double innerProd(const CAssyVector *pX) const;
	void updateCommBoundary();
	void sumupCommBoundary();

	int restrictTo(CAssyVector *pVc) const;
	int prolongateFrom(const CAssyVector *pVc);

	size_t getNumOfVector() const;
	const CVector *getVector(size_t index) const;
	CVector *getVector(size_t index);
	
	void print_elem() const {
		for (std::vector<CVector*>::const_iterator iv = mvVector.begin(); iv != mvVector.end(); iv++) (*iv)->print_elem();
	};

        uint& getDOF();

        void dump();//2011.01.12 列ベクトルのダンプ

private:
	CAssyModel *mpAssyModel;

	std::vector<CVector*> mvVector;
	// commTable

        uint mnDOF;//アセンブル方程式のDOF
};

}

#endif /* ASSYVECTOR_H_ */
