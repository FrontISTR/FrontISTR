/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/MatrixBCRS.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "TypeDef.h"
#include "Logger.h"

#include "Matrix.h"
#include <vector>
#include <map>
using namespace std;

#include <boost/format.hpp> //fdump_mat

#include <boost/tuple/tuple.hpp>
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost;
using namespace boost::numeric;

#ifdef INTELMKL
#include <mkl_boost_ublas_matrix_prod.hpp>
#endif

namespace pmw
{
#ifndef MATRIXBCRS_H_
#define MATRIXBCRS_H_
class CMesh;
class CMatrixBCRS: public CMatrix
{
public:
    CMatrixBCRS(CMesh *pMesh, const uiint& nDOF);
    virtual ~CMatrixBCRS();

    uiint Matrix_Add_Nodal(const uiint& inode, const uiint& jnode, const double* NodalMatrix);
    uiint Matrix_Add_Elem(CMesh *pMesh, const uiint& iElem, double *ElemMatrix);

    void Matrix_Clear();

    virtual void multVector(CVector *pV, CVector *pP) const;

    void setValue_D(const uiint& inode,const uiint& idof, const double& value);
    void setValue(const uiint& inode, const uiint& idof, const double& dDiag, CVector *pRHS, const double& dRHS);
    void setZero_NonDiag(const uiint& inode, const uiint& idof, CVector *pRHS, const double& dRHS);

    uiint setupPreconditioner(iint type);
    uiint setupSmoother(iint type);

    double inverse(ublas::matrix<double> pA, ublas::matrix<double> *pB);
    double determinant(ublas::matrix<double> pA);
    void transpose(ublas::matrix<double> pA, ublas::matrix<double> *pB);
    void print_elem(ublas::matrix<double> pA);

    uiint precond(const CVector *pR, CVector *pZ) const;
    uiint relax(const CVector *pF, CVector *pV) const;
    //--
    // MG : Restriction 行列
    //--
    void restrict(CMatrixBCRS *pCrsMat);


    uiint& getDOF(){ return mnDOF;}

    void dump();


    void initOverlapComm();
    void increment_OverlapComm(const uiint& inode);

protected:
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

    iint mPrecond;

    vector<uint8>  mvNumOverlapComm;//---- 節点毎の通信所属領域オーバーラップ数


    //----------------------
    //  matrix restrict 操作
    //----------------------
    CMesh *mpMesh;// '13.02.18 追加
    
    void VecMarge(vuint& v);//---- vectorマージ
    //--
    // RA :単純CRS サイズ：C×F
    //--
    vuint mvRAIndex,mvRAItem;
    vdouble mvRA;//----- 配列数 = 要素数*mnDOF*mnDOF
    map<uiint, map<uiint, uiint> > mmRCIndex;// [row][col]= index : 行・列 番号からRAの配列Indexを取得

    //--
    // R :単純CRS サイズ：C×F
    //--
    vuint mvRIndex,mvRItem;
    vdouble mvR;
    vdouble mvAggValue;//------リストリクトのAggregate値
    
    // 行列×行列の一部 : 行×列の積 : rowVecは自由度任意、colVecは自由度1
    vdouble Prod_RowCol(vdouble& rowVec, uiint& nDOF, vuint& rowItem, vdouble& colVec, vuint& colItem);
    bool existItem(vuint& vItem, uiint id, uiint& pos);// Item列にIDが存在するか. 存在する場合の配列位置

    //--
    // ij-map : mvAL, mvAU のi行-j列に対応するindex
    //--
    map<uiint, map<uiint,uiint> > mALHash, mAUHash;

    void fdump_R(uiint nNumCrsNode, vuint& vIndex, vuint& vItem, vdouble& vR);
    void fdump_RA(uiint nNumCrsNode, vuint& vIndex, vuint& vItem, vdouble& vRA, uiint& nDOF);
    void fdump_mat(uiint nLevel);
    //
    // matrix restrict 操作 end
    //-------------------------
};
#endif /* MATRIXBCRS_H_ */
}

