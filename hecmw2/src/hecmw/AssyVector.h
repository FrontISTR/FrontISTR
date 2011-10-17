/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/AssyVector.h
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
	const CVector::ElemType &operator[](uiint idx) const;
	CVector::ElemType &operator[](uiint idx);
	const double &operator()(uiint meshID, uiint nodeID, uiint dof) const;
	double &operator()(uiint meshID, uiint nodeID, uiint dof);
        void Vector_Clear(const uiint& iMesh);
	void setZero();
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
        void dump();
private:
	CAssyModel *mpAssyModel;
	std::vector<CVector*> mvVector;
        uiint mnDOF;
};
}
#endif /* ASSYVECTOR_H_ */
