/*
 * MatrixBCRS.h
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
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

class CMatrixBCRS: public CMatrix
{
public:
	CMatrixBCRS(CMesh *pMesh, const uiint& nDOF);
	virtual ~CMatrixBCRS();

        uiint Matrix_Add_Nodal(const uiint& inode, const uiint& jnode, const double* NodalMatrix);
	uiint Matrix_Add_Elem(CMesh *pMesh, const uiint& iElem, double *ElemMatrix);
        
        void Matrix_Clear();// Matrix 0 clear (非線形の行列更新の準備)
        
	virtual void multVector(CVector *pV, CVector *pP) const;

	void setValue_D(const uiint& inode,const uiint& idof, const double& value);//単純に対角項の値をセット

        void setValue(const uiint& inode, const uiint& idof, const double& dDiag, CVector *pRHS, const double& dRHS); //対角項に値をセット、mvD[inode]の行と列を0
        void setZero_NonDiag(const uiint& inode, const uiint& idof, CVector *pRHS, const double& dRHS);

	uiint setupPreconditioner(iint type);
	uiint setupSmoother(iint type);

	double inverse(ublas::matrix<double> pA, ublas::matrix<double> *pB);
	double determinant(ublas::matrix<double> pA);
	void transpose(ublas::matrix<double> pA, ublas::matrix<double> *pB);
	void print_elem(ublas::matrix<double> pA);

	uiint precond(const CVector *pR, CVector *pZ) const;
	uiint relax(const CVector *pF, CVector *pV) const;

        void dump();//2011.01.05 行列をダンプ
private:
	uiint mnNode;
	uiint mnNodeInternal;
	uiint mnDOF;
	uiint mINL;
	uiint mINU;
	std::vector<uiint> mvIndexL;
	std::vector<uiint> mvIndexU;
	std::vector<uiint> mvItemL;
	std::vector<uiint> mvItemU;
	std::vector<ublas::matrix<double> > mvD;
	std::vector<ublas::matrix<double> > mvAL;
	std::vector<ublas::matrix<double> > mvAU;
	std::vector<ublas::matrix<double> > mvALU;
	std::vector<double> mvWW;

////        std::vector<std::vector<double> > mvdAikX;// 0セットした行列の右辺への振替値: 配列数：Node数×DOF数
        
	iint mPrecond;
};

}

#endif /* MATRIXBCRS_H_ */
