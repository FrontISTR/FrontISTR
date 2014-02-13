/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/EqnTerm.h
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
namespace pmw
{
#ifndef EQNTERM_H_
#define EQNTERM_H_
class CEqnTerm
{
public:
    CEqnTerm() {}
    virtual ~CEqnTerm() {}

    //--
    // Term:データ
    //--
    void set(uiint meshID, uiint meshIX, uiint nodeID, uiint nodeIX, uiint dof, double coef, uiint comm_rank, uiint solv_rank, uiint overlap_num, vuint& vRank) {
        mMeshID = meshID;
        mNodeID = nodeID;

        mMeshIX = meshIX;//--- '12.04.06
        mNodeIX = nodeIX;//--- '12.04.06

        mDOF = dof;
        mCoef= coef;

        mCommRank= comm_rank;
        mSolvRank= solv_rank;

        mOverlapNum = overlap_num;
        mvRank= vRank;
    }
    const uiint& meshID() const {
        return mMeshID;
    }
    const uiint& nodeID() const {
        return mNodeID;
    }

    const uiint& meshIX() const {
        return mMeshIX;   //--- '12.04.06
    }
    const uiint& nodeIX() const {
        return mNodeIX;   //--- '12.04.06
    }

    const uiint& dof() const {
        return mDOF;
    }
    const double& coef() const {
        return mCoef;
    }

    const uiint& comm_rank() const {
        return mCommRank;
    }
    const uiint& solv_rank() const {
        return mSolvRank;
    }

    const uiint& overlap_num() const {
        return mOverlapNum;
    }

    const vuint& overlap_rank() const {
        return mvRank;
    }

private:
    uiint mMeshID;
    uiint mNodeID;
    uiint mMeshIX;//--- '12.04.06
    uiint mNodeIX;//--- '12.04.06
    uiint mDOF;
    double mCoef;

    uiint mCommRank;//--通信ランク
    uiint mSolvRank;//--演算ランク

    uiint mOverlapNum;//----- MPC接合面の領域分割オーバーラップ数
    vuint mvRank;//---------- MPC接合面の領域分割オーバーラップ・ランク
};
#endif /* EQNTERM_H_ */
}
