/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/MatrixBCRS.cpp
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
#include "MatrixBCRS.h"
#include "Mesh.h"
#include "Node.h"
#include "Element.h"
#include "Vector.h"
#include <algorithm>
#ifdef _OPENMP
#include "omp.h"
#include "time.h"
#endif
namespace pmw
{
CMatrixBCRS::CMatrixBCRS(CMesh *pMesh, const uiint& nDOF)
{
#ifdef ADVANCESOFT_DEBUG
    printf("enter CMatrixBCRS::CMatrixBCRS \n");
#endif
    mpMesh = pMesh;// '13.02.18

    mnDOF = nDOF;
    mnNode = pMesh->getNumOfNode();
    CIndexBucket *pBucket= pMesh->getBucket();
    mvIndexL.resize(mnNode+1);
    mvIndexU.resize(mnNode+1);
    mvIndexL[0] = 0;
    mvIndexU[0] = 0;
    for (uiint i_node = 0; i_node < mnNode; i_node++) {
        std::vector<uiint> v_item_l;
        std::vector<uiint> v_item_u;
        CNode *pNode= pMesh->getNodeIX(i_node);
        uiint i_node_id = pNode->getID();
        CElement *pElement;
        CAggregateElement *pAggElement= pMesh->getAggElem(i_node_id);
        uiint nNumOfElement= pAggElement->getNumOfElement();

        //cout << "MatrixBCRS construct nNumOfElement:" << nNumOfElement << endl;

        uiint i_elem;
        for(i_elem=0; i_elem < nNumOfElement; i_elem++){
            pElement= pAggElement->get(i_elem);
            uiint k_node, nNumOfNode= pElement->getNumOfNode();

            //cout << "MatrixBCRS construct nNumOfNode:" << nNumOfNode << endl;

            for (k_node = 0; k_node < nNumOfNode; k_node++) {
                uiint k_node_id = pElement->getNode(k_node)->getID();
                uiint k_index = pBucket->getIndexNode(k_node_id);
                if(k_node_id < i_node_id){
                    v_item_l.push_back(k_index);
                }else if(i_node_id < k_node_id){
                    v_item_u.push_back(k_index);
                }
            };
        };
        std::sort(v_item_l.begin(), v_item_l.end());
        std::vector<uiint>::iterator new_end = std::unique(v_item_l.begin(), v_item_l.end());
        v_item_l.erase(new_end, v_item_l.end());
        mvIndexL[i_node + 1] = mvIndexL[i_node] + v_item_l.size();
        mvItemL.insert(mvItemL.end(), v_item_l.begin(), v_item_l.end());

        std::sort(v_item_u.begin(), v_item_u.end());
        new_end = std::unique(v_item_u.begin(), v_item_u.end());
        v_item_u.erase(new_end, v_item_u.end());
        mvIndexU[i_node + 1] = mvIndexU[i_node] + v_item_u.size();
        mvItemU.insert(mvItemU.end(), v_item_u.begin(), v_item_u.end());
    }
    mvD.resize(mnNode);
    for (uiint i = 0; i < mnNode; i++) {
        mvD[i].resize(mnDOF, mnDOF);
        for(uiint i1=0; i1<mnDOF; i1++) for(uiint i2=0; i2<mnDOF; i2++) mvD[i](i1,i2)=0.0;
    }
    mvALU.resize(mnNode);
    for (uiint i = 0; i < mnNode; i++) {
        mvALU[i].resize(mnDOF, mnDOF);
        for(uiint i1=0; i1<mnDOF; i1++) for(uiint i2=0; i2<mnDOF; i2++) mvALU[i](i1,i2)=0.0;
    }
    mINL = mvIndexL[mnNode];
    mvAL.resize(mINL);
    for (uiint i = 0; i < mINL; i++) {
        mvAL[i].resize(mnDOF, mnDOF);
        for(uiint i1=0; i1<mnDOF; i1++) for(uiint i2=0; i2<mnDOF; i2++) mvAL[i](i1,i2)=0.0;
    }
    mINU = mvIndexU[mnNode];
    mvAU.resize(mINU);
    for (uiint i = 0; i < mINU; i++) {
        mvAU[i].resize(mnDOF, mnDOF);
        for(uiint i1=0; i1<mnDOF; i1++) for(uiint i2=0; i2<mnDOF; i2++) mvAU[i](i1,i2)=0.0;
    }

    //--
    // ij-map
    //--
    for(uiint row=0; row < mnNode; row++){
        for(uiint i=mvIndexL[row]; i < mvIndexL[row+1]; i++){
            uiint col=mvItemL[i];
            mALHash[row][col] = i;
        };
        for(uiint i=mvIndexU[row]; i < mvIndexU[row+1]; i++){
            uiint col=mvItemU[i];
            mAUHash[row][col] = i;
        };
    };

    //--
    // overlap number
    //--
    initOverlapComm();//-------------- 初期化:通信オーバーラップ数
    uiint nNumOfCommMesh= pMesh->getCommMesh2Size();
    for(uiint icom=0; icom < nNumOfCommMesh; icom++){
        CCommMesh2 *pCommMesh2= pMesh->getCommMesh2IX(icom);

        uiint nNumOfCommNode= pCommMesh2->getCommNodeSize();
        for(uiint icnode=0; icnode < nNumOfCommNode; icnode++){
            CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);

            uiint node_id= pCommNode->getNodeID();
            uiint inode= pBucket->getIndexNode(node_id);

            increment_OverlapComm(inode);//---- 通信オーバーラップ数++
        };

        //////debug
        ////cout << "MatrixBCRS コンストラクタ icom:" << icom << endl;
        ////for(uiint icnode=0; icnode < nNumOfCommNode; icnode++){
        ////    cout << " icnode:" << icnode << "  overlap_num:" << (int)mvNumOverlapComm[icnode] << endl;
        ////};
    };

#ifdef ADVANCESOFT_DEBUG	
    for (uiint i = 0; i < mnNode; i++) {
        printf(" %ld %ld %ld ; ", i, mvIndexL[i], mvIndexL[i+1]-1);
        for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
                printf(" %ld",mvItemL[j]);
        }
        cout << endl;
    }
    for (uiint i = 0; i < mnNode; i++) {
        printf(" %ld %ld %ld ; ", i, mvIndexU[i], mvIndexU[i+1]-1);
        for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
                printf(" %ld",mvItemU[j]);
        }
        cout << endl;
    }
    printf("exit CMatrixBCRS::CMatrixBCRS \n");
#endif
}
CMatrixBCRS::~CMatrixBCRS()
{
}
uiint CMatrixBCRS::Matrix_Add_Nodal(const uiint& inode, const uiint& jnode, const double* NodalMatrix)
{
    uiint nMatSize = mnDOF;
    uiint kL, kU;
    bool bL,bU;
    uiint irow = inode;
    uiint icol = jnode;
    if( icol == irow ) {
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            mvD[icol](ii,jj) += NodalMatrix[nMatSize*ii+jj];
        };
    } else if( icol < irow ) {
        bL = false;
        for (uiint k = mvIndexL[irow]; k < mvIndexL[irow+1]; k++) {
            if( icol == mvItemL[k] ) {
                kL = k;
                bL = true;
            }
        };
        if( !bL ) printf("***** error in matrix index ***** %ld %ld %s \n",irow,icol,"-1");
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            mvAL[kL](ii,jj) += NodalMatrix[nMatSize*ii+jj];
        };
    } else if( icol > irow ) {
        bU = false;
        for (uiint k = mvIndexU[irow]; k < mvIndexU[irow+1]; k++) {
            if( icol == mvItemU[k] ) {
                kU = k;
                bU = true;
            }
        };
        if( !bU ) printf("***** error in matrix index ***** %ld %ld %s \n",irow,icol,"-1");
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            mvAU[kU](ii,jj) += NodalMatrix[nMatSize*ii+jj];
        };
    }
    return 1;
}
uiint CMatrixBCRS::Matrix_Add_Elem(CMesh *pMesh, const uiint& iElem, double *ElemMatrix)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CMatrixBCRS::Matrix_Add_Elem %ld %e \n", iElem, ElemMatrix[0]);
#endif
    CElement *pElement = pMesh->getElementIX(iElem);
    vector<CNode*> vNode= pElement->getNode();
    uiint nLocalNode = vNode.size();
    uiint nMatSize = nLocalNode * mnDOF;
    CIndexBucket *pBucket = pMesh->getBucket();
    for(uiint i=0; i< nLocalNode; i++) for(uiint j=0; j< nLocalNode; j++){
        uiint kL, kU;
        bool bL, bU;
        uiint iNodeID = vNode[i]->getID();//------------ node id
        uiint jNodeID = vNode[j]->getID();//------------ node id
        uiint irow = pBucket->getIndexNode(iNodeID);//-- node index
        uiint icol = pBucket->getIndexNode(jNodeID);//-- node index
        if( icol == irow ) {
            for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                mvD[icol](ii,jj) += ElemMatrix[nMatSize*(ii+i*mnDOF) + jj+j*mnDOF];
            };
        } else if( icol < irow ) {
            bL = false;
            for (uiint k = mvIndexL[irow]; k < mvIndexL[irow+1]; k++) {
                if( icol == mvItemL[k] ) { 
                    kL = k;
                    bL = true;
                }
            };
            if( !bL ) printf("***** error in matrix index ***** %ld %ld %ld %ld %s \n",i,j,irow,icol,"-1");
            for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                mvAL[kL](ii,jj) += ElemMatrix[nMatSize*(ii+i*mnDOF)+jj+j*mnDOF];

                ////mvAL[kL](ii,jj) /= mvNumOverlapComm[icol];//--------------- オーバーラップ領域数による除算
            };
        } else if( icol > irow ) {
            bU = false;
            for (uiint k = mvIndexU[irow]; k < mvIndexU[irow+1]; k++) {
                if( icol == mvItemU[k] ) { 
                    kU = k;
                    bU = true;
                }
            };
            if( !bU ) printf("***** error in matrix index ***** %ld %ld %ld %ld %s \n",i,j,irow,icol,"-1");
            for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                mvAU[kU](ii,jj) += ElemMatrix[nMatSize*(ii+i*mnDOF) + jj+j*mnDOF];

                        ////mvAU[kU](ii,jj) /= mvNumOverlapComm[icol];//--------------- オーバーラップ領域数による除算
                    }
            };
        }
#ifdef ADVANCESOFT_DEBUG
    printf(" exit CMatrixBCRS::Matrix_Add_Elem \n");
#endif
    return 1;
}
void CMatrixBCRS::Matrix_Clear()
{
    for(uiint i=0; i< mnNode; i++) {
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                mvD[i](ii,jj) = 0.0;
            };
    };
    for(uiint i=0; i < mINL; i++) {
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                mvAL[i](ii,jj) = 0.0;
            };
    };
    for(uiint i=0; i < mINU; i++) {
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                mvAU[i](ii,jj) = 0.0;
            };
    };
}
//--
// A*v=P
//--
void CMatrixBCRS::multVector(CVector *pV, CVector *pP) const
{
#ifdef _OPENMP
    #pragma omp parallel
    for (uiint i = omp_get_thread_num() * mnNode / omp_get_num_threads(); i < (omp_get_thread_num() + 1) * mnNode; i++) {
#else
    for (uiint i = 0; i < mnNode; i++) {
#endif
        (*pP)[i] = prod(mvD[i], (*pV)[i]);
        for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
            (*pP)[i] += prod(mvAL[j], (*pV)[mvItemL[j]]);
        };
        for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
            (*pP)[i] += prod(mvAU[j], (*pV)[mvItemU[j]]);
        };
    };
}
void CMatrixBCRS::setValue_D(const uiint& inode, const uiint& idof, const double& value)
{
    mvD[inode](idof, idof) = value;
}
//--
// Dirichlet 直接代入(2/2)
//--
void CMatrixBCRS::setValue(const uiint& inode, const uiint& idof, const double& dDiag, CVector *pRHS, const double& dRHS)
{
    for(uiint j=0; j < mnDOF; j++) {
        if(idof != j) {
            mvD[inode](idof, j) = 0.0;//------------ inode-idof 行を"0.0"で払う
        }
    };
    for(uiint j=0; j < mnDOF; j++) {
        if(idof != j) {
            double diffVal = -mvD[inode](j, idof) * dRHS;//-- 減算
            pRHS->addValue(inode, j, diffVal);//--------- 右辺を減算(0で払った行列値*境界値)
            mvD[inode](j, idof) = 0.0;//----------------- mvD の非対角項のidof列を"0.0"
        }
    };
    mvD[inode](idof,idof) = dDiag;//---------------- 対角項を1.0
}
//--
// Dirichlet 直接代入(1/2)
//--
void CMatrixBCRS::setZero_NonDiag(const uiint& inode, const uiint& idof, CVector* pRHS, const double& dRHS)
{
    //--
    // 行を0で払い、右辺は境界値 : 対角項=1で解は境界値になる.
    //--
    for(uiint j=mvIndexL[inode]; j < mvIndexL[inode+1]; j++) {
        for(uiint jdof=0; jdof < mnDOF; jdof++) {
            mvAL[j](idof,jdof) = 0.0;//-------------- 行は”0.0”で払うだけ.
        };
    };
    for(uiint j=mvIndexU[inode]; j < mvIndexU[inode+1]; j++) {
        for(uiint jdof=0; jdof < mnDOF; jdof++) {
            mvAU[j](idof,jdof) = 0.0;//-------------- 行は”0.0”で払うだけ.
        };
    };

    //--
    // 列を0で払い、右辺から行列要素*境界値を減算 : 左辺との整合性をとる.
    //--
    for(uiint jnode=0; jnode < mnNode; jnode++) {
        for(uiint i=mvIndexL[jnode]; i < mvIndexL[jnode+1]; i++) {
            uiint col=mvItemL[i];
            //列がinode
            if(col==inode) {
                for(uiint jdof=0; jdof < mnDOF; jdof++) {
                    double diffVal = -mvAL[i](jdof,idof) * dRHS;//--- 減算
                    pRHS->addValue(jnode, jdof, diffVal);
                    mvAL[i](jdof,idof) = 0.0;//-------------------------- 列を0.0で払う
                };
            }
        };
    };
    for(uiint jnode=0; jnode < mnNode; jnode++) {
        for(uiint i=mvIndexU[jnode]; i < mvIndexU[jnode+1]; i++) {
            uiint col=mvItemU[i];
            //列がinode
            if(col==inode) {
                for(uiint jdof=0; jdof < mnDOF; jdof++) {
                    double diffVal = -mvAU[i](jdof, idof) * dRHS;//--- 減算
                    pRHS->addValue(jnode, jdof, diffVal);
                    mvAU[i](jdof,idof) = 0.0;//-------------------------- 列を0.0で払う
                };
            }
        };
    };
}
//--
//
//--
void CMatrixBCRS::VecMarge(vuint& v)
{
    sort(v.begin(), v.end());
    vuint::iterator new_end = unique(v.begin(), v.end());
    v.erase(new_end, v.end());
}
//--
// MG:制限
//--
void CMatrixBCRS::restrict(CMatrixBCRS* pCrsMat)
{
    if(pCrsMat==NULL) return;//----- Level:0 (コースグリッド)

    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    //--
    // # コースグリッド・ノード 周囲への操作なので、情報を集め直す
    //   TODO: => コンストラクタへ移動予定
    //--
    uiint nLevel= mpMesh->getMGLevel();//カレントMG_Level
    uiint nNumOfCGrid= pCrsMat->mnNode;

    CIndexBucket *pBucket= mpMesh->getBucket();

    // R (サイズ：C×F) : 単純CRS  
    mvRIndex.resize(nNumOfCGrid+1);
    mvRIndex[0]=0;

    //--
    // R : コンストラクタ => A行列のコースグリッド節点周囲(係数)と同じ
    //--
    for(uiint inode=0; inode < nNumOfCGrid; inode++){//---コースグリッド・ノード Loop
        CNode* pNode= mpMesh->getNodeIX(inode);
        vuint vENodeIX;//----------pNode周囲のNodeIX:FEMコネクティビティ
        
        if(pNode->getNumOfParentNode(nLevel) != 0)//----------Error
            pLogger->Info(Utility::LoggerMode::Error, " MatrixBCRS::restrict, NumOfParentNode != 0");

        uiint nCoreNodeIX = pBucket->getIndexNode(pNode->getID());//-----節点IX(コースグリッド)
        vENodeIX.push_back(nCoreNodeIX);

        // 1. Aggregate Value 計算
        // 2. restrict係数の取得
        CAggregateElement* pAggElem= mpMesh->getAggElemIX(inode);
        //
        // 一旦、周囲の節点IX番号を集める.
        //
        for(uiint ielem=0; ielem < pAggElem->getNumOfElement(); ielem++){
            CElement *pElem= pAggElem->get(ielem);

            for(uiint enode=0; enode < pElem->getNumOfNode(); enode++){
                CNode* pENode= pElem->getNode(enode);

                // pNodeを除外 : 既に取得済み
                if(pENode->getID() != pNode->getID() ){
                    uiint nNodeIX= pBucket->getIndexNode(pENode->getID());//----節点IX(周辺グリッド)
                    vENodeIX.push_back(nNodeIX);
                }

            };//for(elem node)end
        };//for(agg elem)end


        VecMarge( vENodeIX );//------- vectorマージ


        double AggVal=0.0;
        vdouble vCoeff;
        // pNode周囲のNodeを集める
        // pNode周囲のAggValueを求める
        for(uiint around=0; around < vENodeIX.size(); around++){
            uiint ix= vENodeIX[around];
            CNode* pArndNode= mpMesh->getNodeIX(ix);

            uiint nPareNum= pArndNode->getNumOfParentNode(nLevel);

            if(nPareNum!=0){
                double Coeff = 1.0/(double)nPareNum;
                AggVal += Coeff;
                vCoeff.push_back(Coeff);
            }else{
                //--
                // pNode自身 : Coeff=1.0
                //--
                double Coeff = 1.0;
                AggVal += Coeff;
                vCoeff.push_back(Coeff);
            }
        };//for(around)end

        //--
        // R : 係数を利用する時にはAggValueで除算 => 最初から除算するとPで利用できない.
        //--
        mvAggValue.push_back(AggVal);//----------- pNode周囲の係数加算値

        mvRItem.insert(mvRItem.end(), vENodeIX.begin(), vENodeIX.end());//------ pNode周囲のNodeIX(列番号)
        mvRIndex[inode + 1]= mvRIndex[inode] + vENodeIX.size();//----------行の先頭位置保存
        mvR.insert(mvR.end(), vCoeff.begin(), vCoeff.end());//-------------R : pNode周囲の係数(自由度は無関係)
        
    };//for( CoarseGridノード ) end


    ////cout << " MatrixBCRS::restrict ---- R constructed" << endl;

    ////
    ////fdump_R(nNumOfCGrid, mvRIndex, mvRItem, mvR);
    ////

    //--
    // RA 確保(サイズ：C×F) : コンストラクタ
    // # サイズは R と同じ(DOF=1の場合)
    // # mmRCIndex : Hash-Hash [row][col]= index : 行・列 番号にからRAのIndexを取得
    //--
    mvRAIndex= mvRIndex;//------ R と同じIndexのはず
    mvRAItem= mvRItem;//-------- R と同じItem のはず
    uiint nMaxI = mvRIndex[nNumOfCGrid];//--MaxIndex
    mvRA.resize(nMaxI*mnDOF*mnDOF);
    // RA 初期化
    for(uiint i=0; i < nMaxI*mnDOF*mnDOF; i++){
        mvRA[i]=0.0;
    };
    // mmRCIndex 初期化
    for(uiint inode=0; inode < nNumOfCGrid; inode++){
        for(uiint i=mvRAIndex[inode]; i < mvRAIndex[inode+1]; i++){
            uiint col= mvRAItem[i];
            mmRCIndex[inode][col]= i;//----------- 行・列 番号に対応するRAの配列インデックス
        };
    };
    
    ///////////// debug
    ////ofstream ofs;
    ////ofs.open("dump_calc_RA.txt");
    /////////////

    //--
    // RA = R * A  : Rの行とAの列を掛ける => Rの列番号:Aの行番号
    //--
    // # inode 周囲のRestrict
    //--
    for(uiint inode=0; inode < nNumOfCGrid; inode++){
        for(uiint i=mvRIndex[inode]; i < mvRIndex[inode+1]; i++){
            uiint col_r= mvRItem[i];//-----------------------------これに対応する、Aの行 : Rの列番号がAの行番号
            
            ////ofs << " inode:" << inode << " col_r:" << col_r << endl;
            
            //
            // 下三角
            //
            for(uiint pos_a=mvIndexL[col_r]; pos_a < mvIndexL[col_r+1]; pos_a++){
                uiint col_a = mvItemL[pos_a];//------------------- col_r行のAの列番号

                ////ofs << " AL 列:" << col_a << endl;

                for(uiint idof=0; idof < mnDOF; idof++){
                for(uiint jdof=0; jdof < mnDOF; jdof++){
                    //
                    // # Rのcol列に対応する、Aのcol行
                    //
                    uiint i_ra = mmRCIndex[inode][col_a];//------------------------- RAのインデックス(R行番号・A列番号)
                    uiint nPos = (i_ra * mnDOF*mnDOF) + (idof*mnDOF) + jdof;//----- 自由度を含めたRAのindex番号

                    mvRA[nPos] += mvR[i]*mvAL[pos_a](idof,jdof)/mvAggValue[inode];//-- RA = R * A
                };
                };
            };//col_r行の列Item

            //
            // 上三角
            //
            for(uiint pos_a=mvIndexU[col_r]; pos_a < mvIndexU[col_r+1]; pos_a++){
                uiint col_a = mvItemU[pos_a];//-------------------------- col_r行のAの列番号

                ////ofs << " AU 列:" << col_a << endl;

                for(uiint idof=0; idof < mnDOF; idof++){
                for(uiint jdof=0; jdof < mnDOF; jdof++){
                    //
                    // # Rのcol列 => Aのcol行
                    //
                    uiint i_ra = mmRCIndex[inode][col_a];//------------------------ RAのインデックス(R行番号・A列番号)
                    uiint nPos = (i_ra * mnDOF*mnDOF) + (idof*mnDOF) + jdof;//----- 自由度を含めたRAのindex番号

                    mvRA[nPos] += mvR[i]*mvAU[pos_a](idof,jdof)/mvAggValue[inode];//-- RA = R * A
                };
                };
            };//col_r行の列Item

            //
            // R の[inode行][col列] の index => i
            //
            for(uiint idof=0; idof < mnDOF; idof++){
            for(uiint jdof=0; jdof < mnDOF; jdof++){
                // 対角
                // # Rの列番号:Aの行番号
                uiint i_ra = mmRCIndex[inode][col_r];
                uiint nPos = (i_ra * mnDOF*mnDOF) + (idof*mnDOF) + jdof;

                mvRA[nPos] += mvR[i]*mvD[col_r](idof,jdof)/mvAggValue[inode];//-- RA = R * A
            };
            };
            
        };//for(index):Rのインデックス
    };//for(inode) 行:Rの行

    ////////////////////
    ////ofs.close();// debug
    ////////////////////

    ////cout << " MatrixBCRS::restrict ---- RA constructed" << endl;

    /////////////////////////////////////////////////////
    ////fdump_RA(nNumOfCGrid, mvRAIndex, mvRAItem, mvRA, mnDOF);
    ////fdump_mat(nLevel);//-- ファイングリッドMat
    //////////////////////////////////////////////////////

    //--
    // コースグリッドAc : ゼロ・クリア
    //--
    pCrsMat->Matrix_Clear();


    ////cout << " MatrixBCRS::restrict ---- pCrsMat cleard" << endl;
    
    ////////////////////////////
    ////ofs.open("dump_calc_Ac.txt");
    ////////////////////////////

    //--
    //  Ac = RA * P  :  P(サイズ：F×C) => 単純CCSは、CRSの転置なので”R"をそのまま利用する.
    //--
    for(uiint inode=0; inode < nNumOfCGrid; inode++){
        uiint is=mvRAIndex[inode], ie=mvRAIndex[inode+1];
        uiint nRowSize= (ie-is)*mnDOF*mnDOF;
        vdouble vRow;
        vRow.resize( nRowSize );
        //
        // vRow <=: mvRAから一行代入
        //
        uiint nDOF2 = mnDOF*mnDOF;
        for(uiint i=is; i < ie; i++){
            uiint ir= i-is;
            for(uiint idof=0; idof < nDOF2; idof++){
                vRow[ir*nDOF2 + idof] = mvRA[i*nDOF2 + idof];
            };
        };
        vuint vRowItem;
        vRowItem.resize(ie-is);
        for(uiint i=is; i < ie; i++){
            uiint ir= i-is;
            vRowItem[ir] = mvRAItem[i];
        };
        //
        // P は、R 転置なので、Rの行をPの列として扱う
        //    =>  CRSとして作ったものをCCS扱いする
        //
        for(uiint jnode=0; jnode < nNumOfCGrid; jnode++){
            uiint js=mvRIndex[jnode], je=mvRIndex[jnode+1];
            vdouble vCol;
            uiint nColSize = je-js;
            vCol.resize( nColSize );
            //
            // vCol <=: vRから一行代入
            //
            for(uiint i=js; i < je; i++){
                uiint ic = i-js;
                vCol[ic] = mvR[i];
            };
            vuint vColItem;
            vColItem.resize(nColSize);
            for(uiint i=js; i < je; i++){
                uiint ic= i-js;
                vColItem[ic] = mvRItem[i];
            };
            //--
            // 1行 × 1列 積
            //--
            vdouble vVal= Prod_RowCol( vRow, mnDOF, vRowItem, vCol, vColItem);
            
            //--
            // Ac[inode行][jnode列] = mvRa[i]*vR[j];
            //--
            uiint index;
            // 下三角
            if(inode > jnode){
                index = mALHash[inode][jnode];
                for(uiint idof=0; idof < mnDOF; idof++){
                for(uiint jdof=0; jdof < mnDOF; jdof++){
                    uiint nPos = idof*mnDOF + jdof;
                    pCrsMat->mvAL[index](idof,jdof) = vVal[nPos];
                };
                };
            }
            // 対角
            if(inode==jnode){
                for(uiint idof=0; idof < mnDOF; idof++){
                for(uiint jdof=0; jdof < mnDOF; jdof++){
                    uiint nPos = idof*mnDOF + jdof;
                    pCrsMat->mvD[inode](idof,jdof) = vVal[nPos];
                };
                };
            }
            // 上三角
            if(inode < jnode){
                index = mAUHash[inode][jnode];
                for(uiint idof=0; idof < mnDOF; idof++){
                for(uiint jdof=0; jdof < mnDOF; jdof++){
                    uiint nPos = idof*mnDOF + jdof;
                    pCrsMat->mvAU[index](idof,jdof) = vVal[nPos];
                };
                };
            }
        };//for(jnode):Pの列

    };//for(inode) 行:RAの行

    //////////////
    ////ofs.close();
    //////////////

    ////cout << " MatrixBCRS::restrict ---- pCrsMat constructed" << endl;

    //////////////////////////////
    ////pCrsMat->fdump_mat(nLevel-1);
    //////////////////////////////
}
//--
// 行列×行列の一部 : 行×列の積 : rowVecは自由度任意、colVecは自由度1
//--
vdouble CMatrixBCRS::Prod_RowCol(vdouble& rowVec, uiint& nDOF, vuint& rowItem, vdouble& colVec, vuint& colItem)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    
    ////cout << " MatrixBCRS::Prod_RowCol, rowVec.size:" << rowVec.size() << " colVec.size:" << colVec.size() << endl;
    ////cout << " MatrixBCRS::Prod_RowCol, rowItem.size:" << rowItem.size() << " colItem.size:" << colItem.size() << endl;

    //--
    // vVal : DOF×DOF行列を一列にしたベクトル
    //--
    //  0 1 2
    //  3 4 5  => 0 1 2 3 4 5 6 7 8
    //  6 7 8

    vdouble vVal; vVal.resize(nDOF*nDOF);
    for(uiint i=0; i < nDOF*nDOF; i++) vVal[i]=0.0;

    ////cout << " MatrixBCRS::Prod_RowCol, rowItem.size:" << rowItem.size() << endl;

    //--
    // 行 進行
    //--
    for(uiint i=0; i < rowItem.size(); i++){
        uiint rowID= rowItem[i];//行Vecの{ 値の存在する列番号 }
        uiint posCol;

        ////cout << " MatrixBCRS::Prod_RowCol, i:" << i << endl;

        //
        // 行番号に対応する列番号が存在した場合に積をとる:なければ”0”
        //
        if( existItem(colItem, rowID, posCol) ){

            ////cout << " MatrixBCRS::Prod_RowCol, rowID:" << rowID << " posCol:" << posCol << endl;

            for(uiint idof=0; idof < nDOF; idof++){
            for(uiint jdof=0; jdof < nDOF; jdof++){

                uiint posRow = (i*nDOF*nDOF) + (idof*nDOF) + jdof;
                uiint posVal = (idof*nDOF) + jdof;

                vVal[posVal] += rowVec[posRow]*colVec[posCol];
            };
            };
        }// if (existItem) end
    };
    
    return vVal;
}
//--
// Item列にIDが存在するか. 存在する場合の配列位置
//--
bool CMatrixBCRS::existItem(vuint& vItem, uiint id, uiint& pos)
{
    uiint low=0;
    uiint high=vItem.size()-1;

    uiint ix;
    while( low <= high ){
        ix=(low+high)/2;

        if( id == vItem[ix] ){
            pos=ix;
            return true;
        }else if( id < vItem[ix] ){
            if(ix!=0){
                high=ix-1;
            }else{
                return false;//---- 負は無い.
            }
        }else{
            low=ix+1;
        }
    };
    return false;
}
//--
// R のファイル出力(単純CRSの出力)
//--
void CMatrixBCRS::fdump_R(uiint nNumCrsNode, vuint& vIndex, vuint& vItem, vdouble& vR)
{
    ofstream ofs;
    ofs.open("drmp_R.txt");

    ofs << "---- R ----" << endl;
    for(uiint inode=0; inode < nNumCrsNode; inode++){
        for(uiint i=vIndex[inode]; i < vIndex[inode+1]; i++){
            uiint ix= vItem[i];
            double val= vR[i];
            //////////////////////////////////////////////////////////////////
            ofs << " row:" << inode << " col:" << ix << " val:" << val << endl;
            //////////////////////////////////////////////////////////////////
        };
    };
    ofs.close();
}
//--
// RA のファイル出力(単純CRSの出力) : DOF*DOF
//--
void CMatrixBCRS::fdump_RA(uiint nNumCrsNode, vuint& vIndex, vuint& vItem, vdouble& vRA, uiint& nDOF)
{
    ofstream ofs;
    ofs.open("dump_RA.txt");

    ofs << "---- RA ----" << endl;
    for(uiint inode=0; inode < nNumCrsNode; inode++){
        for(uiint i=vIndex[inode]; i < vIndex[inode+1]; i++){
            uiint ix= vItem[i];
            ofs << " row:" << inode << " col:" << ix << endl;
            
            for(uiint idof=0; idof < nDOF; idof++){
                for(uiint jdof=0; jdof < nDOF; jdof++){
                    uiint nPos= (i*nDOF*nDOF) + (idof*nDOF) + jdof;
                    double val= vRA[nPos];
                    ////////////////////////////////
                    ofs << format(" %1$15.8e") %val ;
                    ////////////////////////////////
                };
                ofs << endl;
            };
        };
    };
    ofs.close();
}
//--
// matrix ファイル出力
//--
void CMatrixBCRS::fdump_mat(uiint nLevel)
{
    string sLevel;
    stringstream ss;
    ss << nLevel;
    ss >> sLevel;

    string sfnameD("dump_D_"+sLevel+".txt");
    string sfnameL("dump_L_"+sLevel+".txt");
    string sfnameU("dump_U_"+sLevel+".txt");

    ofstream ofsD,ofsL,ofsU;
    ofsD.open(sfnameD.c_str());
    ofsL.open(sfnameL.c_str());
    ofsU.open(sfnameU.c_str());

    ofsD << " ---- mvD ---- " << endl;
    for(uiint i=0; i< mnNode; i++){
        ofsD << "row:" << i << " col:" << i << endl;
        for(uiint ii=0; ii < mnDOF; ii++){
            for(uiint jj=0; jj < mnDOF; jj++) {
                ///////////////////////////////////////////
                ofsD << format(" %1$15.8e")  %mvD[i](ii,jj);
                ///////////////////////////////////////////
            };
            ofsD << endl;
        };
    };
    ofsD.close();

    ofsL << " ---- mvAL ---- " << endl;
    for(uiint inode=0; inode < mnNode; inode++){
        for(uiint i=mvIndexL[inode]; i < mvIndexL[inode+1]; i++){
            uiint ix=mvItemL[i];
            ofsL << "row:" << inode << " col:" << ix << endl;
            for(uiint ii=0; ii < mnDOF; ii++){
                for(uiint jj=0; jj < mnDOF; jj++) {
                    ////////////////////////////////////////////
                    ofsL << format(" %1$15.8e")  %mvAL[i](ii,jj);
                    ////////////////////////////////////////////
                };
                ofsL << endl;
            };
        };
    };
    ofsL.close();

    ofsU << " ---- mvAU ---- " << endl;
    for(uiint inode=0; inode < mnNode; inode++){
        for(uiint i=mvIndexU[inode]; i < mvIndexU[inode+1]; i++){
            uiint ix=mvItemU[i];
            ofsU << "row:" << inode << " col:" << ix << endl;
            for(uiint ii=0; ii < mnDOF; ii++){
                for(uiint jj=0; jj < mnDOF; jj++) {
                    ///////////////////////////////////////////
                    ofsU << format(" %1$15.8e") %mvAU[i](ii,jj);
                    ///////////////////////////////////////////
                };
                ofsU << endl;
            };
        };
    };
    ofsU.close();

}



//--
//
//--
uiint CMatrixBCRS::setupPreconditioner(iint type)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CMatrixBCRS::setupPreconditioner \n");
#endif
    ublas::matrix<double> pA(mnDOF,mnDOF), pB(mnDOF,mnDOF), pC(mnDOF,mnDOF);
    mPrecond = type;
    iint itype;
    itype = -99;
    if( mPrecond == 1 ) itype = 2;
    if( mPrecond == 2 ) itype = 2;
    if( mPrecond == 3 ) itype = 2;
    if( mPrecond == 4 ) itype = 1;
    if( itype == -99 ) printf("setup precondition; fatal error %ld\n", mPrecond);
    switch( itype ) {
    case( 1 ):
        printf(" [TYPE:1] Preconditioner \n");
        for (int i = 0; i < mnNode; i++) {
            pA = mvD[i];
            for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
                transpose(mvAL[j], &pB);
                pC  = prod(mvAL[j], mvALU[mvItemL[j]] );
                pA -= prod(pC, pB);
            };
            inverse(pA, &mvALU[i]);
        };
        break;
    case( 4 ):
        printf(" [TYPE:4] Preconditioner \n");
        for (uiint i = 0; i < mnNode; i++) mvALU[i] = mvD[i];
        for (uiint i = 0; i < mnNode; i++) {
            pA = mvALU[i];
            inverse(pA, &mvALU[i]);
            pA = mvALU[i];
            for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
                pC  = prod(pA, mvAU[j]);
                transpose(mvAU[j], &pB);
                mvALU[mvItemU[j]] -= prod(pB, pC);
            };
        };
        break;
    case( 2 ):
        printf(" [TYPE:2] Preconditioner \n");
        for (uiint i = 0; i < mnNode; i++) {
            inverse(mvD[i], &mvALU[i]);
        };
        break;
    case( 3 ):
        printf(" [TYPE:3] Preconditioner \n");
        break;
    }
#ifdef ADVANCESOFT_DEBUG
    printf(" exit CMatrixBCRS::setupPreconditioner \n");
#endif
    return 1;
}
double CMatrixBCRS::inverse(ublas::matrix<double> pA, ublas::matrix<double> *pB)
{
    double det = determinant( pA );
    double Recip = 1.0 / det;
    (*pB)(0, 0) = Recip * ( pA(1, 1) * pA(2, 2) - pA(2, 1) * pA(1, 2) );
    (*pB)(0, 1) = Recip * (-pA(0, 1) * pA(2, 2) + pA(2, 1) * pA(0, 2) );
    (*pB)(0, 2) = Recip * ( pA(0, 1) * pA(1, 2) - pA(1, 1) * pA(0, 2) );
    (*pB)(1, 0) = Recip * (-pA(1, 0) * pA(2, 2) + pA(2, 0) * pA(1, 2) );
    (*pB)(1, 1) = Recip * ( pA(0, 0) * pA(2, 2) - pA(2, 0) * pA(0, 2) );
    (*pB)(1, 2) = Recip * (-pA(0, 0) * pA(1, 2) + pA(1, 0) * pA(0, 2) );
    (*pB)(2, 0) = Recip * ( pA(1, 0) * pA(2, 1) - pA(2, 0) * pA(1, 1) );
    (*pB)(2, 1) = Recip * (-pA(0, 0) * pA(2, 1) + pA(2, 0) * pA(0, 1) );
    (*pB)(2, 2) = Recip * ( pA(0, 0) * pA(1, 1) - pA(1, 0) * pA(0, 1) );
    return det;
}
void CMatrixBCRS::transpose(ublas::matrix<double> pA, ublas::matrix<double> *pB)
{
    for (int i=0; i< mnDOF; i++) for (int j=0; j<mnDOF; j++) (*pB)(i, j) = pA(j, i);
}
double CMatrixBCRS::determinant(ublas::matrix<double> pA)
{
    double det = pA(0, 0) * pA(1, 1) * pA(2, 2)
                 + pA(1, 0) * pA(2, 1) * pA(0, 2)
                 + pA(2, 0) * pA(0, 1) * pA(1, 2)
                 - pA(2, 0) * pA(1, 1) * pA(0, 2)
                 - pA(1, 0) * pA(0, 1) * pA(2, 2)
                 - pA(0, 0) * pA(2, 1) * pA(1, 2);
    return det;
}
uiint CMatrixBCRS::setupSmoother(iint type)
{
    printf(" enter CMatrixBCRS::setupSmoother \n");
    printf(" exit CMatrixBCRS::setupSmoother \n");
    return 1;
}

//--
// 前処理 & 前処理MGのsmoother
//--
uiint CMatrixBCRS::precond(const CVector *pR, CVector *pZ) const
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    CVector::ElemType WW(mnDOF);//------------ ublas::vector
    ublas::zero_vector<double> vzero(mnDOF);

    iint itype;
    itype = -99;
    if( mPrecond == 1 ) itype = 1;//1:Jacobi
    if( mPrecond == 2 ) itype = 2;//2:MG, MG前処理の場合は、SSORを内部の前処理とする.
    //  if( mPrecond == 2 ) itype = 1;//test 11.08.01
    if( mPrecond == 3 ) itype = 2;//3:SSOR
    if( mPrecond == 4 ) itype = 1;//4:ILU, ILU前処理の場合は、Jacobiとする.

    //---- 1〜4以外の前処理(mPrecond)が選択された場合Error
    if( itype == -99 ) printf("precond; fatal error %ld\n", mPrecond, itype);

    switch( itype ) {
    case( 0 ):
        pZ->subst(pR);
        break;
    case( 1 )://---- Jacobi
        pZ->subst(pR);
        for(uiint i=0; i< mnNode; i++) {
            for(uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++)  (*pZ)[i] -= prod( mvAL[j], (*pZ)[mvItemL[j]] );

            (*pZ)[i] = prod(mvALU[i], (*pZ)[i]);

////////////            // "nan" 判定
////////////            double val;
////////////            for(uiint idof=0; idof < mnDOF; idof++){
////////////                val= (*pZ)[i](idof);
////////////                if(isnan(val))
////////////                    pLogger->Info_format(Utility::LoggerMode::Error, "%s%d%s%d", "MatrixBCRS::precond, inode:", i, " idof:", idof);
////////////            };
        };
        for(uiint i= mnNode-1; i >= 0 && i < mnNode; i--) {

            for(uiint idof=0; idof < mnDOF; idof++)   WW(idof)=0.0;
            for(uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++)   WW += prod( mvAU[j], (*pZ)[mvItemU[j]] );

            (*pZ)[i] -= prod(mvALU[i], WW);
        };
        break;
    case( 2 )://---- SSOR
        for(uiint loop=0; loop < 10; loop++) {
            #pragma omp parallel for private(WW)
            for(uiint i=0; i< mnNode; i++) {
                WW = (*pR)[i];
                for(uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
                    WW -= prod( mvAL[j], (*pZ)[mvItemL[j]] );
                };
                for(uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
                    WW -= prod( mvAU[j], (*pZ)[mvItemU[j]] );
                };
                WW = prod( mvALU[i], WW );// WW = invD * WW

                double omega = 0.5;//-------係数
                (*pZ)[i] = (*pZ)[i] + omega * ( WW - (*pZ)[i] );
            };// i-mnNode end
        };// loop-10 end
        break;
    case( 3 )://---- relax:SSOR
        for(uiint i=0; i< mnNode; i++) (*pZ)[i] = vzero;
        for(uiint i=0; i< 50; i++) relax(pR, pZ);
        break;
    }
    return 1;
}
//
// MG_Solver本体のsmoother : pR:右、pZ:解
//
uiint CMatrixBCRS::relax(const CVector *pR, CVector *pZ) const
{
    CVector::ElemType WW(mnDOF);

    //---- SSOR
    for(uiint i=0; i< mnNode; i++) {
        WW = (*pR)[i];
        for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
            WW -= prod(mvAL[j], (*pZ)[mvItemL[j]]);
        };
        for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
            WW -= prod(mvAU[j], (*pZ)[mvItemU[j]]);
        };
        WW = prod(mvALU[i], WW);
        double omega = 1.8;//--------係数(加速)
        (*pZ)[i] = (*pZ)[i] + omega * ( WW - (*pZ)[i] );
    };

////    //---- Jacobi
////    pZ->subst(pR);
////    for(uiint i=0; i< mnNode; i++) {
////        for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
////            (*pZ)[i] -= prod(mvAL[j], (*pZ)[mvItemL[j]]);
////        };
////        (*pZ)[i] = prod(mvALU[i], (*pZ)[i]);
////    };
////    for(uiint i= mnNode-1; i >= 0 && i < mnNode; i--){
////
////        for(uiint idof=0; idof < mnDOF; idof++) WW(idof)=0.0;
////
////        for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
////            WW += prod(mvAU[j], (*pZ)[mvItemU[j]]);
////        };
////        (*pZ)[i] -= prod(mvALU[i], WW);
////    };

    return 1;
}


void CMatrixBCRS::dump()
{
    cout << " ---- mvD ---- " << endl;
    for(uiint i=0; i< mnNode; i++) {
        cout << "Node i:" << i << " ";
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                cout << mvD[i](ii,jj) << "  ";
            };
        cout << endl;
    };
    cout << " ---- mvAL ---- " << endl;
    for(uiint i=0; i < mINL; i++) {
        cout << "mINL:" << i << " ";
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                cout << mvAL[i](ii,jj) << "  ";
            };
        cout << endl;
    };
    cout << " ---- mvAU ---- " << endl;
    for(uiint i=0; i < mINU; i++) {
        cout << "mINU:" << i << " ";
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                cout << mvAU[i](ii,jj) << "  ";
            };
        cout << endl;
    };
}
void CMatrixBCRS::initOverlapComm()
{
    mvNumOverlapComm.resize(mnNode);

    for(uiint inode=0; inode < mnNode; inode++) {
        mvNumOverlapComm[inode]= 1;//------------- 最低1つの領域には属している.
    };
}
void CMatrixBCRS::increment_OverlapComm(const uiint& inode)
{
    mvNumOverlapComm[inode]++;
}

}//namespace pmw


