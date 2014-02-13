/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/AssyVector.h
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
#include "Logger.h"
#include <vector>
#include "Vector.h"
#include "HEC_MPI.h"
namespace pmw
{
#ifndef ASSYVECTOR_H_
#define ASSYVECTOR_H_
class CAssyModel;

class CAssyVector
{
public:
    CAssyVector();

    void initAssyVector(CAssyModel *pAssyModel, vuint& vDOF);//--vDOF:グローバルMesh番号(MeshID)に対応するDOF
    void initAssyVector(const CAssyVector *pAssyVector);

    CAssyVector(CAssyModel *pAssyModel, vuint& vDOF);//----------vDOF:グローバルMesh番号(MeshID)に対応するDOF
    CAssyVector(const CAssyVector *pAssyVector);

    virtual ~CAssyVector();

    uiint size() const;
    uiint size();

    uiint getNumOfMesh() const;
    uiint getNumOfMesh();

    uiint& getDOF(const uiint& imesh) const;
    uiint& getDOF(const uiint& imesh);

    const CVector::ElemType &operator[](uiint idx) const;
    CVector::ElemType &operator[](uiint idx);
    const double &operator()(uiint meshID, uiint nodeID, uiint dof) const;
    double &operator()(uiint meshID, uiint nodeID, uiint dof);
    void Vector_Clear(const uiint& iMesh);
    void Vector_Clear();
    void setZero();
    void setValue(uiint imesh, uiint inode, uiint idof, double value);
    void addValue(const uiint& imesh, const uiint& inode, const uiint& idof, const double& value);
    double& getValue(uiint imesh, uiint inode, uiint idof);

    void set(const CAssyVector *pV);
    void add(const CAssyVector *pV);
    void subst(const CAssyVector *pV);
    double norm2() const;
    double innerProd(const CAssyVector *pX) const;

    void subtrac(const CAssyVector *pV);//減算

    void updateCommBoundary();//-- 通信(1/2)
    void sumupCommBoundary(); //-- 通信(2/2)


    uiint restrictTo(CAssyVector *pVc) const;    //-- MGリストリクション
    uiint prolongateFrom(const CAssyVector *pVc);//-- MGプロロンゲーション

    uiint getNumOfVector() const;
    const CVector *getVector(uiint index) const;
    CVector *getVector(uiint index);

    void print_elem() const {
        for (std::vector<CVector*>::const_iterator iv = mvVector.begin(); iv != mvVector.end(); iv++) (*iv)->print_elem();
    };

    void dump() const;
    void dump();

    //--
    // MGCycle テスト関数
    //--
    void setCoarseVector(CAssyVector *coarse) {
        mpCoarseVector = coarse;
    };
    CAssyVector* getCoarseVector() {
        return mpCoarseVector;
    };
    CAssyVector* getCoarseVector() const {
        return mpCoarseVector;
    };

    void divis_OverlapNum();//---- ベクトル要素を通信オーバーラップ数で除算

private:
    CAssyModel *mpAssyModel;
    std::vector<CVector*> mvVector;

    CAssyVector *mpCoarseVector;//------ コースグリッドAssyVector
};
#endif /* ASSYVECTOR_H_ */
}

