/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Vector.cpp
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
#include "Vector.h"
#include "Mesh.h"
#include "Node.h"
namespace pmw
{
CVector::CVector(CMesh *pMesh, const uiint& nDOF)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CVector::CVector \n");
#endif
    mpMesh = pMesh;
    mnDOF = nDOF;
    mnNode = pMesh->getNumOfNode();
    mvVector.resize(mnNode);
    for (uiint i = 0; i < mnNode; i++) {
        mvVector[i].resize(mnDOF);
        // 初期化: 0.0
        for(uiint idof=0; idof < mnDOF; idof++) mvVector[i](idof)=0.0;
    }

    initOverlapComm();//--- 通信所属領域数の初期化

#ifdef ADVANCESOFT_DEBUG
    printf(" exit CVector::CVector \n");
#endif
}
// コピーコンストラクター
CVector::CVector(const CVector *pVector)
{
    mnDOF = pVector->mnDOF;
    mnNode = pVector->mnNode;
    mnNodeInternal = pVector->mnNodeInternal;
    mvVector.resize(mnNode);

    for(uiint inode = 0; inode < mnNode; inode++) {
        mvVector[inode].resize(mnDOF);
    };

    //--
    // 通信オーバーラップ数の初期化 & 代入
    //--
    initOverlapComm();

    for(uiint inode=0; inode < mnNode; inode++) {
        mvNumOverlapComm[inode]= pVector->getOverlapComm_Num(inode);
    };
}
CVector::~CVector()
{
    //--mvVectorは実体なのでDeleteせず.
}
uiint CVector::size() const
{
    return mnNode;
}
bool CVector::isScopeNode(const uiint& idx) const
{
    uiint max = mvVector.size()-1;
    if(idx > max) {
        return false;
    } else {
        return true;
    }
}
const CVector::ElemType& CVector::operator[](uiint idx) const
{
    return mvVector[idx];
}
CVector::ElemType& CVector::operator[](uiint idx)
{
    return mvVector[idx];
}
void CVector::Vector_Clear()
{
    uiint i,idof;
    for(i=0; i < mnNode; i++) {
        for(idof=0; idof < mnDOF; idof++) {
            //mvVector[i][idof] = 0.0;
            mvVector[i](idof) = 0.0;
        };
    };
}
void CVector::setZero()
{
    for (uiint i = 0; i < mnNode; i++) {
        mvVector[i].clear();
    }
}
void CVector::setValue(uiint inode, uiint idof, double value)
{
    //cout << "CVector::setValue  inode:" << inode << " idof:" << idof << " value:" << value << endl;

    //mvVector[inode][idof] = value;
    mvVector[inode](idof) = value;
}
void CVector::addValue(uiint inode, uiint idof, double value)
{
    //mvVector[inode][idof] += value;
    mvVector[inode](idof) += value;
}
double& CVector::getValue(uiint inode, uiint idof)
{
    //return( mvVector[inode][idof] );
    return mvVector[inode](idof);
}
void CVector::sumSV(double alpha, const CVector *pX, CVector *pY) const
{
    for (uiint i=0; i < mnNode; i++) {
        pY->mvVector[i] = mvVector[i] + alpha * pX->mvVector[i];
    }
}
void CVector::addSV(double alpha, const CVector *pX)
{
    for (uiint i=0; i < mnNode; i++) {
        mvVector[i] += alpha * pX->mvVector[i];
    }
}

// 代入
void CVector::set(const CVector* pX)
{
    for (uiint i=0; i < mnNode; i++) {
        mvVector[i] = pX->mvVector[i];
    }
}
// 加算
void CVector::add(const CVector *pX)
{
    for (uiint i=0; i < mnNode; i++) {
        mvVector[i] += pX->mvVector[i];
    }
}
// 代入
void CVector::subst(const CVector *pX)
{
    for (uiint i=0; i < mnNode; i++) {
        mvVector[i] = pX->mvVector[i];
    }
}
double CVector::norm2() const
{
    return innerProd(this);
}
double CVector::innerProd(const CVector *pX) const
{
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)

    for (uiint i = 0; i < mnNode; i++) {
        for (uiint j = 0; j < mnDOF; j++) {
            sum +=  mvVector[i](j) * pX->mvVector[i](j);
        }
    }
    return sum;
}

//減算
void CVector::subtrac(const CVector *pX)
{
    for(uiint i=0; i < mnNode; i++) {
        mvVector[i] -= pX->mvVector[i];
    }
}

void CVector::updateCommBoundary()
{
}


void CVector::VecMarge(vuint& v) const
{
    sort(v.begin(), v.end());
    vuint::iterator new_end = unique(v.begin(), v.end());
    v.erase(new_end, v.end());
}

//
// 残差制限(Restrict)
//
uiint CVector::restrictTo(CVector *pV) const // pV:コースグリッド
{
    uiint nLevel= mpMesh->getMGLevel();

    CIndexBucket *pBucket= mpMesh->getBucket();
    CAggregateElement *pAggElem;
    uiint nCGridSize= pV->size();//---------- コースグリッドNode数
    //
    // コースグリッド節点へ制限(restrict)
    //
    for( uiint i=0; i< nCGridSize; i++) {
        
        if(mpMesh->isDirichletBNode(i)) continue;//--- Dirichlet境界は除外
        
        vuint vArndNode;
        ////CNode* pNode = mpMesh->getNodeIX(i);
        pAggElem= mpMesh->getAggElemIX(i);

        uiint nNumE = pAggElem->getNumOfElement();
        // 周囲要素
        for(uiint ielem=0; ielem < nNumE; ielem++){
            CElement *pElem = pAggElem->get(ielem);

            uiint nNumArndNode = pElem->getNumOfNode();
            for(uiint enode=0; enode < nNumArndNode; enode++){
                CNode *pArndNode = pElem->getNode(enode);
                uiint ix = pBucket->getIndexNode(pArndNode->getID());
                vArndNode.push_back(ix);
            };
        };
        //--
        // vector マージ
        //--
        VecMarge(vArndNode);

        //--
        // restrict value.
        //--
        double dCoeff(0.0);// coeff集合値
        vdouble vVal;
        vVal.resize(mnDOF);
        for(uiint idof=0; idof < mnDOF; idof++) vVal[idof]=0.0;

        for(uiint iarnd=0; iarnd < vArndNode.size(); iarnd++){
            uiint ix= vArndNode[iarnd];
            //
            // Dirichlet境界以外のNode値を集める:ファイングリッドで境界になっているNodeを除外
            //
            if( !mpMesh->isDirichletBNode(ix) ){
                CNode *pArndNode = mpMesh->getNodeIX(ix);
                uiint nNumPare = pArndNode->getNumOfParentNode(nLevel);

                if(nNumPare!=0){
                    dCoeff += 1.0/nNumPare;
                    for(uiint idof=0; idof < mnDOF; idof++)
                        vVal[idof] += (1.0/nNumPare)*mvVector[ix](idof);
                }else{
                    dCoeff += 1.0;
                    for(uiint idof=0; idof < mnDOF; idof++)
                        vVal[idof] += mvVector[ix](idof);// i==ix のはず.
                }
            }//---- if(!isBNode)end
        };

        for(uiint idof=0; idof < mnDOF; idof++)
            (*pV)[i](idof) = vVal[idof]/dCoeff;//---- pV:コースグリッド

    };// for(コースグリッド) end

    return 0;
}
//
// 修正量伸展(Prolongate)
//
uiint CVector::prolongateFrom(const CVector *pV) // pV:コースグリッド
{
    CIndexBucket *pBucket= mpMesh->getBucket();
    uiint nLevel= mpMesh->getMGLevel();

    for( uiint i=0; i< mpMesh->getNumOfNode(); i++) {
        CNode* pNode = mpMesh->getNodeIX(i);
        uiint nNumP = pNode->getNumOfParentNode(nLevel);

        // 加算するので、一旦"0"クリア
        for(uiint idof=0; idof < mnDOF; idof++) mvVector[i](idof) = 0.0;//--- 0 clear
        
        if( nNumP==0 ) {
            //
            // 子Node==親Node : pV(コースグリッド)
            //
            //----
            // コースグリッドがディレクレ境界の場合は加算しない.
            //----
            if( !mpMesh->isDirichletBNode(i) )
                for(uiint idof=0; idof < mnDOF; idof++) mvVector[i](idof)= (*pV)[i](idof);
            
        } else {
            //
            // 子Node==複数の親Node : pV(コースグリッド)
            //
            for(uiint j=0; j < nNumP; j++) {
                CNode* pPareNode= pNode->getParentNode(nLevel,j);
                uiint kID = pPareNode->getID();
                uiint k = pBucket->getIndexNode(kID);
                //----
                //コースグリッドがディレクレ境界の場合は加算しない.
                // ! 連続番号 なので、ファイングリッドでも、コースグリッド数内であれば、コースグリッドの境界判定可能
                //----
                if( !mpMesh->isDirichletBNode(k) )
                    for(uiint idof=0; idof < mnDOF; idof++) mvVector[i](idof) += (*pV)[k](idof)/(double)nNumP;
            };
        }//if(nNumP==0)

    };// for(ファイングリッド)end

    return 0;
}




void CVector::dump()
{
    uiint i,j;

////    // DOF:0 のみ出力
////    for(i = 0; i < mnNode; i++){
////        cout << mvVector[i](0) << " ";
////    };
////    cout << endl;

    // 全DOFを出力
    for(j=0; j < mnDOF; j++) {
        cout << " DOF:" << j << " : ";
        for(i=0; i < mnNode; i++) {
            cout << mvVector[i](j) << " ";
        };
        cout << endl;
    };
}
//--
// 節点毎の通信所属領域オーバーラップ '12.03.23
//--
void CVector::initOverlapComm()
{
    mvNumOverlapComm.resize(mnNode);

    for(uiint inode=0; inode < mnNode; inode++) {
        mvNumOverlapComm[inode]= 1;//------------- 最低1つの領域には属している.
    };
}
void CVector::increment_OverlapComm(const uiint& inode)
{
    mvNumOverlapComm[inode]++;
}
uint8 CVector::getOverlapComm_Num(const uiint& inode) const
{
    return mvNumOverlapComm[inode];
}
//--
// ベクトル要素をオーバーラップ数で除算
//--
void CVector::divis_OverlapNum()
{
    ////CHecMPI *pMPI= CHecMPI::Instance();
    ////cout << " Vector::divis_OverlapNum  ---- enter  rank:" << pMPI->getRank() << endl;

    for(uiint inode=0; inode < mnNode; inode++) {
        double dDiv = 1.0/mvNumOverlapComm[inode];

        ////cout << " Vector::divis_OverlapNum inode:" << inode << " dDiv:"<< dDiv << " rank:" << pMPI->getRank() << endl;//debug

        for(uiint idof=0; idof < mnDOF; idof++) {
            mvVector[inode][idof] *= dDiv;
        };
    };

    ////cout << " Vector::divis_OverlapNum  ---- exit   rank:" << pMPI->getRank() << endl;
}

}// namespace END

