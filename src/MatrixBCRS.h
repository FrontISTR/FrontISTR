/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   MatrixBCRS.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef MATRIXBCRS_H_
#define MATRIXBCRS_H_
#include "TypeDef.h"
#include "Matrix.h"
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric;
namespace pmw
{
class CMesh;
class CSolverPre;
class CMatrixBCRS: public CMatrix
{
public:
	CMatrixBCRS(/* const */ CMesh *pMesh);
	virtual ~CMatrixBCRS();
	int Matrix_Add_Elem(CMesh *pMesh, uint iElem, double *ElemMatrix);
	virtual void multVector(CVector *pV, CVector *pP) const;
	void multVector111(const CVector *pV, CVector *pP) const;
	void setValue(int inode, int idof, double value);
	void MatrixMult33(ublas::matrix<double> A, std::vector<double> B, std::vector<double> C);
	int setupPreconditioner(int type);
	int setupSmoother(int type);
	double inverse(ublas::matrix<double> pA, ublas::matrix<double> *pB);
	double determinant(ublas::matrix<double> pA);
	void transpose(ublas::matrix<double> pA, ublas::matrix<double> *pB);
	void print_elem(ublas::matrix<double> pA);
	int precond(const CVector *pR, CVector *pZ) const;
	int relax(const CVector *pF, CVector *pV) const;
private:
	int mnNode;
	int mnNodeInternal;
	int mnDOF;
	int mINL;
	int mINU;
	std::vector<int> mvIndexL;
	std::vector<int> mvIndexU;
	std::vector<int> mvItemL;
	std::vector<int> mvItemU;
	std::vector<ublas::matrix<double> > mvD;
	std::vector<ublas::matrix<double> > mvAL;
	std::vector<ublas::matrix<double> > mvAU;
	std::vector<ublas::matrix<double> > mvALU;
	std::vector<double> mvWW;
	int mPrecond;
	friend class CSolverPreILU;
};
}
#endif /* MATRIXBCRS_H_ */
