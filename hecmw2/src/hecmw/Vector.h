/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Vector.h
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
#ifndef VECTOR_H_
#define VECTOR_H_
#include "TypeDef.h"
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#ifdef INTELMKL
#include <mkl_boost_ublas_matrix_prod.hpp>
#endif
using namespace boost::numeric;

#define SOL_VEC 0
#define RHS_VEC 1

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

    void set(const CVector *pX);//---- 代入
    void add(const CVector *pX);//---- 加算
    void subst(const CVector *pX);//-- 代入
    double norm2() const;
    double innerProd(const CVector *pX) const;

    void subtrac(const CVector *pX);//-- 減算

    void updateCommBoundary();

    uiint restrictTo(CVector *pV) const;
    uiint prolongateFrom(const CVector *pV);

    void print_elem() const ;
    void dump();

    uiint& getDOF() {
        return mnDOF;
    }

    //-- 節点毎の通信所属領域オーバーラップ '12.03.23
    void initOverlapComm();
    void increment_OverlapComm(const uiint& inode);//------ 通信オーバーラップ数の++
    uint8 getOverlapComm_Num(const uiint& inode) const;//-- 通信オーバーラップ数:CVectorのコピーコンストラクタ用

    void divis_OverlapNum();//--- ベクトル要素をオーバーラップ数で除算
private:
    uiint mnNode;
    uiint mnNodeInternal;
    uiint mnDOF;
    std::vector<ElemType> mvVector;//[inode][idof]:  ElemTypeは, ublas_vector: DOF別の値を格納
    CMesh *mpMesh;
    bool isScopeNode(const uiint& idx) const;

    vector<uint8>  mvNumOverlapComm;//---- 節点毎の通信所属領域オーバーラップ数

        void VecMarge(vuint& v ) const;//---- vectorマージ
};
}
#endif /* VECTOR_H_ */
