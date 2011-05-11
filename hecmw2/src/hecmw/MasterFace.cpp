
#include "ContactNode.h"

//
//  MasterFace.cpp
//
//  MPCマスター面群の面
//			2009.10.15
//			2009.01.08
//			k.Takeda
#include <vector>

#include "Vertex.h"

#include "MasterFace.h"
using namespace pmw;


// construct & destruct
//
CMasterFace::CMasterFace()
{
    mvVector.resize(3);//内外判定時のスレーブ点-頂点2点で構成される面の面ベクトルの一時的な入れ物
    
    // L1とL2ベクトルの最近接点を表すパラメータa,b,c,d,e,f, 最近接点を表すスカラー値s,t
    // --
    mS=0.0; mT=0.0;
    mParam_A=0.0; mParam_B=0.0; mParam_C=0.0; mParam_D=0.0; mParam_E=0.0; mParam_F=0.0;
    mvParam_R.resize(3); mvParam_R[0]=0.0; mvParam_R[1]=0.0; mvParam_R[2]=0.0;

    //    mNzVector.resize(3);//正規化ベクトル計算時の入れ物
}
CMasterFace::~CMasterFace()
{
    //std::cout << "~CMasterFace" << std::endl;
}

// CSkinFace::refine()で使用
//
CSkinFace* CMasterFace::generateFace()
{
    mpOtherFace = new CMasterFace;

    return mpOtherFace;
}


// 外積
//   引数:ConNodeを中心とした,頂点2点
//
vdouble& CMasterFace::CrossProduct(CContactNode *pConNode, const uiint& ivert, const uiint& jvert)
{
    CContactNode *pVert0,*pVert1;
    pVert0= mvConNode[ivert];  pVert1= mvConNode[jvert];
    
    // 外積に使用する2つのベクトル
    double X1= pVert0->getX() - pConNode->getX();
    double Y1= pVert0->getY() - pConNode->getY();
    double Z1= pVert0->getZ() - pConNode->getZ();

    double X2= pVert1->getX() - pConNode->getX();
    double Y2= pVert1->getY() - pConNode->getY();
    double Z2= pVert1->getZ() - pConNode->getZ();

    // 外積
    // x = y1 z2 - z1 y2
    // y = z1 x2 - x1 z2
    // z = x1 y2 - y1 x2
    //
    mvVector[0]= Y1*Z2 - Z1*Y2;
    mvVector[1]= Z1*X2 - X1*Z2;
    mvVector[2]= X1*Y2 - Y1*X2;

    return mvVector;
}


// ベクトルの絶対値
//
double CMasterFace::VectorLength(const vdouble& vec)
{
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}
double CMasterFace::VectorLength(CContactNode* pConNode0, CContactNode* pConNode1)
{
    double dX= pConNode1->getX() - pConNode0->getX();
    double dY= pConNode1->getY() - pConNode0->getY();
    double dZ= pConNode1->getZ() - pConNode0->getZ();
    
    return sqrt(dX*dX + dY*dY + dZ*dZ);
}
double CMasterFace::VectorLength(CContactNode* pConNode, const vdouble& vec)
{
    double dX= vec[0]-pConNode->getX();
    double dY= vec[1]-pConNode->getY();
    double dZ= vec[2]-pConNode->getZ();

    return sqrt(dX*dX + dY*dY + dZ*dZ);
}
double CMasterFace::VectorLength(const vdouble& vBasePoint, const vdouble& vPoint)
{
    double dX= vPoint[0]-vBasePoint[0];
    double dY= vPoint[1]-vBasePoint[1];
    double dZ= vPoint[2]-vBasePoint[2];
    
    return sqrt(dX*dX + dY*dY + dZ*dZ);
}
double CMasterFace::VectorLength(const vdouble& vec, CContactNode* pConNode)
{
    double dX= pConNode->getX() - vec[0];
    double dY= pConNode->getY() - vec[1];
    double dZ= pConNode->getZ() - vec[2];

    return sqrt(dX*dX + dY*dY + dZ*dZ);
}

// 線形補間 :ConNode0-ConNode1のidof番目の値の線形補間
//
double CMasterFace::LinearInter(const double& coeff, const uiint& idof,
                                CContactNode* pConBaseNode, CContactNode* pConNode, const uiint& valType)
{
    double Value;
    
    switch(valType){
        case(MPCValueType::Displacement):
            Value= coeff*(pConNode->getDisp(idof) - pConBaseNode->getDisp(idof)) + pConBaseNode->getDisp(idof);
            break;
        case(MPCValueType::Scalar):
            Value= coeff*(pConNode->getScalar(idof) - pConBaseNode->getScalar(idof)) + pConBaseNode->getScalar(idof);
            break;
    }
    return Value;
}
double CMasterFace::LinearInter(const double& coeff, const uiint& idof, const vdouble& vBaseVal, const vdouble& vVal)
{
    return coeff*(vVal[idof]-vBaseVal[idof]) + vBaseVal[idof];
}

// 内外判定のための,ベクトル内積(mvNormalVectorとmvVectorの内積),
//   面ベクトルは,ContactMesh::setupSPointOnMFace()で計算がコール済み.
// --
// mvVectorは,スレーブ点と頂点2点で構成される面の面ベクトル. addSlaveNodeのルーチン内で使用
// --
double& CMasterFace::DotProduct()
{
    mDotValue= mvNormalVector[0]*mvVector[0] + mvNormalVector[1]*mvVector[1] + mvNormalVector[2]*mvVector[2];

    return mDotValue;
}

// 1.与えられたConNodeの"面"内外判定を行って,内ならば追加.
//    ----
//    対象となる点は, 入力データ時点で平面と平面のカップリングが選択してあるのが前提
//      .&&. BBoxにより絞り込みが行われているとする.
//    ----
// 2.スレーブ点が追加される度に,Coefの確保を行う.
//
void CMasterFace::addSlaveNode(CContactNode* pConNode)
{
    int iprodp,iprodm;//外積の+,-のカウンター :: 面ベクトル(正規化)と同じ向きなら”+”. 逆向きなら"-"
    int nNumOfVert;
    uiint ivert,jvert;//頂点Index

    nNumOfVert= this->getNumOfVert();//1次・2次要素どちらでも面の頂点数

    iprodp=0;
    iprodm=0;
    for(ivert=0; ivert< nNumOfVert; ivert++){
        //隣の頂点番号
        if(ivert==nNumOfVert-1){ jvert=0; }else{ jvert= ivert+1; }

        //"pConNodeと頂点2個"の面ベクトル :: mvVector
        CrossProduct(pConNode, ivert, jvert);

        //マスター面の面ベクトルとの内積
        // => 面と同じ方向ならば＋プラス,逆向きならば-マイナス
        if(DotProduct() > 0) iprodp++;
        if(DotProduct() < 0) iprodm--;
        if(DotProduct()== 0){
            iprodp++;//線の上に載っている.'10.02.28
            iprodm--;
        }
    };

    //debug
    //cout << "MasterFace::addSlaveNode, numOfVert= " << nNumOfVert << ", +内外カウント = " << iprodp << ", -内外カウント = " << iprodm << endl;

    //外積のカウンターが(頂点数 or -頂点数)と一致ならば,面内と判定 => Slave点
    if(iprodp==nNumOfVert || iprodm== -nNumOfVert){

        // スレーブ点として追加
        mvSlaveConNode.push_back(pConNode);
        mmSlaveID2Index[pConNode->getID()]= mvSlaveConNode.size()-1;// ID to Index map
        // スレーブ点にマスター面IDを設定
        pConNode->markingSlave();
        pConNode->setMasterFaceID(mID, mLevel);
        
        // Coefの確保
        //  スレーブ点が追加される度にCoefの配列数を増加
        vdouble vCoef;
        uiint nNumOfFaceNode = mvConNode.size();
        vCoef.resize(nNumOfFaceNode);//面のノード数確保(1次・2次 両方の要素に対応するため面のノード数)

        mvvCoef.push_back(vCoef);

    }
    
}

/*
    // L1とL2ベクトルの最近接点を表すパラメータ計算, 最近接点を表すスカラー値s,t
    // --
    double mS,mT, mParam_A, mParam_B, mParam_C, mParam_D, mParam_E, mParam_F;

    double& s_NearCrossP();//L1ベクトルの基点からの距離 s
    double& t_NearCrossP();//L2ベクトルの基点からの距離 t
    // 線形代数方程式の各パラメータ
    void NearCrossP(const vdouble& vL1, const vdouble& vL2, const vdouble& L1P0, const vdouble& L2P0);//L1とL2の最近接点を表す行列要素a,b,c,d,e,fの計算
    // 各パラメータの個別計算関数
    double& a_NearCrossP(const vdouble& vL1);// パラメータa
    double& b_NearCrossP(const vdouble& vL1, const vdouble& vL2);// パラメータb
    double& c_NearCrossP(const vdouble& vL1);// パラメータc
    double& d_NearCrossP();// パラメータd
    double& e_NearCrossP(const vdouble& vL2);// パラメータe
    double& f_NearCrossP(const vdouble& vL2);// パラメータf
    double& r_NearCrossP(const vdouble& L1P0, const vdouble& L2P0);// パラメータr
 */

// NearCrossP(...)でパラメータを求めた後に,最近接点を計算
//
double& CMasterFace::s_NearCrossP()
{
    mS= ((mParam_B * mParam_F)-(mParam_C * mParam_E))/mParam_D;

    return mS;
}
double& CMasterFace::t_NearCrossP()
{
    mT= ((mParam_A * mParam_F)-(mParam_B * mParam_C))/mParam_D;
    
    return mT;
}

// 最近接点のパラメータの一括計算
//
bool CMasterFace::NearCrossP(const vdouble& vL1, const vdouble& vL2, const vdouble& L1P0, CContactNode* pL2P0)
{
    bool bCrossP(true);

    r_NearCrossP(L1P0, pL2P0);

    a_NearCrossP(vL1);
    b_NearCrossP(vL1,vL2);
    c_NearCrossP(vL1);
    d_NearCrossP();
    e_NearCrossP(vL2);
    f_NearCrossP(vL2);

    //debug
    ///////////////////////////////////////
    //cout << "(v,v,v,p) mParam_D= " << mParam_D << endl;
    ///////////////////////////////////////

    if( mParam_D <= 1.0E-6) bCrossP=false;//mParam_Dが0.0ならば平行 => 近接点は存在しない.

    return bCrossP;
}
bool CMasterFace::NearCrossP(const vdouble& vL1, const vdouble& vL2, CContactNode* pL1P0, CContactNode* pL2P0)
{
    bool bCrossP(true);

    r_NearCrossP(pL1P0, pL2P0);

    a_NearCrossP(vL1);
    b_NearCrossP(vL1,vL2);
    c_NearCrossP(vL1);
    d_NearCrossP();
    e_NearCrossP(vL2);
    f_NearCrossP(vL2);

    //debug
    ///////////////////////////////////////
    //cout << "(v,v,p,p) mParam_D= " << mParam_D << endl;
    ///////////////////////////////////////

    if( mParam_D <= 1.0E-6) bCrossP=false;//mParam_Dが0.0ならば平行 => 近接点は存在しない.

    return bCrossP;
}

//// Paramのクリア
//void CMasterFace::clearParam()
//{
//    mParam_A=0.0;
//    mParam_B=0.0;
//    mParam_C=0.0;
//    mParam_D=0.0;
//    mParam_E=0.0;
//    mParam_F=0.0;
//
//    mvParam_R[0]=0.0;
//    mvParam_R[1]=0.0;
//    mvParam_R[2]=0.0;
//}

// NearCrossP()関数の個別パラメータ計算
//
double& CMasterFace::a_NearCrossP(const vdouble& vL1)
{
    return mParam_A= vL1[0]*vL1[0] + vL1[1]*vL1[1] + vL1[2]*vL1[2];
}
double& CMasterFace::b_NearCrossP(const vdouble& vL1, const vdouble& vL2)
{
    return mParam_B= vL1[0]*vL2[0] + vL1[1]*vL2[1] + vL1[2]*vL2[2];
}
double& CMasterFace::c_NearCrossP(const vdouble& vL1)
{
    return mParam_C= vL1[0]*mvParam_R[0] + vL1[1]*mvParam_R[1] + vL1[2]*mvParam_R[2];
}
double& CMasterFace::d_NearCrossP()
{
    return mParam_D= (mParam_A * mParam_E)-(mParam_B * mParam_B);
}
double& CMasterFace::e_NearCrossP(const vdouble& vL2)
{
    return mParam_E= vL2[0]*vL2[0] + vL2[1]*vL2[1] + vL2[2]*vL2[2];
}
double& CMasterFace::f_NearCrossP(const vdouble& vL2)
{
    return mParam_F= vL2[0]*mvParam_R[0] + vL2[1]*mvParam_R[1] + vL2[2]*mvParam_R[2];
}
vdouble& CMasterFace::r_NearCrossP(const vdouble& L1P0, const vdouble& L2P0)
{
    mvParam_R[0]= L1P0[0]-L2P0[0];
    mvParam_R[1]= L1P0[1]-L2P0[1];
    mvParam_R[2]= L1P0[2]-L2P0[2];

    return mvParam_R;
}
vdouble& CMasterFace::r_NearCrossP(CContactNode* pL1P0, CContactNode* pL2P0)
{
    mvParam_R[0]= pL1P0->getX() - pL2P0->getX();
    mvParam_R[1]= pL1P0->getY() - pL2P0->getY();
    mvParam_R[2]= pL1P0->getZ() - pL2P0->getZ();

    return mvParam_R;
}
vdouble& CMasterFace::r_NearCrossP(const vdouble& L1P0, CContactNode* pL2P0)
{
    mvParam_R[0]= L1P0[0] - pL2P0->getX();
    mvParam_R[1]= L1P0[1] - pL2P0->getY();
    mvParam_R[2]= L1P0[2] - pL2P0->getZ();

    return mvParam_R;
}

// 平行ベクトルの場合に,点PへのベクトルをFベクトル辺に投影した距離s (内積:DotProduct)
//
double& CMasterFace::s_ProjecP(const vdouble& vV, CContactNode* pP, CContactNode* pP3)
{
    double dX,dY,dZ;
    dX= pP->getX() - pP3->getX();
    dY= pP->getY() - pP3->getY();
    dZ= pP->getZ() - pP3->getZ();

    mT= vV[0]*dX + vV[1]*dY + vV[2]*dZ;

    return mS;
}
// 平行ベクトルの場合に,点PへのベクトルをEベクトル辺に投影した距離t (内積:DotProduct)
//
double& CMasterFace::t_ProjecP(const vdouble& vV, CContactNode* pP, CContactNode* pP2)
{
    double dX,dY,dZ;
    dX= pP->getX() - pP2->getX();
    dY= pP->getY() - pP2->getY();
    dZ= pP->getZ() - pP2->getZ();
    
    mT= vV[0]*dX + vV[1]*dY + vV[2]*dZ;
    
    return mT;
}

// ベクトルの正規化
void CMasterFace::Normalized(vdouble& vV)
{
    //与えられたベクトルを正規化
    double dLength= sqrt(vV[0]*vV[0] + vV[1]*vV[1] + vV[2]*vV[2]);

    vV[0]/= dLength;
    vV[1]/= dLength;
    vV[2]/= dLength;
}



// EQATION のCoef計算ルーチン
//  Coef = TermA * TermB
//
double CMasterFace::CoefTermA(CContactNode* pOpposNode, const vdouble& inP, CContactNode* pEdgeNode0, CContactNode* pEdgeNode1)
{
    double nume;//numerate 分子
    double deno;//denominator 分母
    double termA;

    nume= VectorLength(pOpposNode, inP);
    deno= VectorLength(pEdgeNode0,pEdgeNode1);

    termA= nume/deno;

    return termA;
}
double CMasterFace::CoefTermB(const vdouble& inP, CContactNode* pSlaveP, const vdouble& inP0, const vdouble& inP1)
{
    double nume;//numerate 分子
    double deno;//denominator 分母
    double termB;

    nume= VectorLength(pSlaveP, inP);
    deno= VectorLength(inP0, inP1);

    termB= nume/deno;

    return termB;
}

// EQUATION  変位
//  ・マスター面の変位 => スレーブ点の変位の計算
//  ・変位は,ContactNodeが所有.
//
// EQUATION  スカラー
//  ・マスター面のスカラー => スレーブ点のスカラーの計算
//  ・スカラー値は,ContactNodeが所有.
//
void CMasterFace::CalcSlave(const uiint& islave, const uiint& valType)
{
    // 
    // * 2直線の最近接点の中間点を近似交点として,マスター面の変位後のスレーブ点を求める.
    //
    
    CContactNode* pSlaveP= mvSlaveConNode[islave];//計算対象のスレーブ点
    
    //// mmvCoef所有の為のID
    //uint slaveID= pSlaveP->getID();

    uiint numOfVert= this->getNumOfVert();

    // 対向する辺の最近接点
    //   E,Fベクトルの最近接点(vNP0,vNP1) => 近似交点(vOP)
    // --
    vdouble vOP;       vOP.resize(3);// E,F ベクトルの近似交点
    vdouble vNP0,vNP1; vNP0.resize(3);vNP1.resize(3);//最近接点2点
    vdouble vE,vF;     vE.resize(3); vF.resize(3);   //E,Fベクトル

    vdouble vPOP;  vPOP.resize(3); //スレーブ点とOP点とのベクトル
    vdouble vP2P3; vP2P3.resize(3);//p2-p3ベクトル
    vdouble vP1P0; vP1P0.resize(3);//p1-p0ベクトル
    vdouble vInterP0,vInterP1; vInterP0.resize(3);vInterP1.resize(3);//p2-p3の補間点, p0-p1の補間点
    double dLenP2P3,dLenP2InP,dLenP0P1,dLenP1InP;//p2-p3の距離,p2-InterP0の距離,p0-p1の距離,p1-InterP1の距離
    vdouble vValInterP0,vValInterP1;//p2-p3補間点の変位値,p0-p1補間点の変位値
    uiint numOfDOF,idof;//ConNodeの自由度数
    
    double dLenInP0InP1,dLenInP0SlaveP;//InterP0-InterP1の距離,InterP0-SlavePの距離
    vdouble vValSlaveP;//スレーブ点での補間値
    
    // 三角形用途
    double  aArea;//三角形全体の面積
    vdouble vArea;//スレーブ点で3分割されてできる,3個の三角形
    vArea.resize(3);

    // パラメータ種類別
    switch(valType){
        case(MPCValueType::Displacement):
            numOfDOF=3;//x,y,zの並進3自由度に固定
            vValInterP0.resize(numOfDOF); vValInterP1.resize(numOfDOF);
            vValSlaveP.resize(numOfDOF);
            break;
        case(MPCValueType::Scalar):
            numOfDOF=mvConNode[0]->getNumOfScalar();
            vValInterP0.resize(numOfDOF); vValInterP1.resize(numOfDOF);
            vValSlaveP.resize(numOfDOF);
            break;
    }

    bool bEF;//左右の辺：EベクトルとFベクトルに交点があるか? の変数

    //debug 用途 ΣCoefのチェック用
    uiint ideb;
    double deb;

    switch(numOfVert){
        case(4)://四辺形
            // E ベクトル p0-p3
            vE[0]= mvConNode[0]->getX() - mvConNode[3]->getX();
            vE[1]= mvConNode[0]->getY() - mvConNode[3]->getY();
            vE[2]= mvConNode[0]->getZ() - mvConNode[3]->getZ();
            // F ベクトル p1-p2
            vF[0]= mvConNode[1]->getX() - mvConNode[2]->getX();
            vF[1]= mvConNode[1]->getY() - mvConNode[2]->getY();
            vF[2]= mvConNode[1]->getZ() - mvConNode[2]->getZ();

            // 正規化
            Normalized(vE);  Normalized(vF);

            // E,F の交点計算 => vOP :E,Fベクトルの最近接点(vOPの計算)
            bEF= NearCrossP(vE,vF, mvConNode[3],mvConNode[2]);
            
            if(bEF){// EとFに最近接点が存在
                s_NearCrossP();// NearCrossPの実行後は,パラメータが揃っているので引数なし
                t_NearCrossP();// NearCrossPの実行後は,パラメータが揃っているので引数なし

                //E,Fベクトルの最近接点2点 NP0,NP1
                vNP0[0]= mS*vE[0] + mvConNode[3]->getX();//ベクトルE * s + ベクトル基点
                vNP0[1]= mS*vE[1] + mvConNode[3]->getY();
                vNP0[2]= mS*vE[2] + mvConNode[3]->getZ();

                vNP1[0]= mT*vF[0] + mvConNode[2]->getX();//ベクトルF * t + ベクトル基点
                vNP1[1]= mT*vF[1] + mvConNode[2]->getY();
                vNP1[2]= mT*vF[2] + mvConNode[2]->getZ();
                
                //E,Fベクトルの近似交点
                vOP[0]= (vNP0[0]+vNP1[0])*0.5;
                vOP[1]= (vNP0[1]+vNP1[1])*0.5;
                vOP[2]= (vNP0[2]+vNP1[2])*0.5;

                ////debug
                ///////////////////////////////////////////
                //{
                //    cout << "スレーブ点 ConID= " << pSlaveP->getID();
                //    cout << "X= " << pSlaveP->getX() << ",Y= " << pSlaveP->getY() << ",Z= " << pSlaveP->getZ() << endl;
                //
                //    cout << " E,Fベクトルの近似交点:" <<
                //            "X= " << vOP[0] << ",Y= " << vOP[1] << ",Z= " << vOP[2] << endl;
                //}
                ///////////////////////////////////////////

                //OPとスレーブ点のベクトル
                vPOP[0]= pSlaveP->getX() - vOP[0];
                vPOP[1]= pSlaveP->getY() - vOP[1];
                vPOP[2]= pSlaveP->getZ() - vOP[2];
                
            }else{// EとFが平行の場合
                ////debug
                //cout << "E-F 平行" << endl;
                
                vPOP[0]= vE[0];// vPOPは,vE,vFと同じ
                vPOP[1]= vE[1];
                vPOP[2]= vE[2];
                
                vOP[0]= pSlaveP->getX();// vOPは,スレーブ点自身
                vOP[1]= pSlaveP->getY();
                vOP[2]= pSlaveP->getZ();
                
            }//if(bEF)の終端
            
            // --
            // InterP0,InterP1の計算, 補間値vValInterP0,vValInterP1の計算
            // --
            
            //p2-p3ベクトル
            vP2P3[0]= mvConNode[3]->getX() - mvConNode[2]->getX();
            vP2P3[1]= mvConNode[3]->getY() - mvConNode[2]->getY();
            vP2P3[2]= mvConNode[3]->getZ() - mvConNode[2]->getZ();
            //p1-p0ベクトル
            vP1P0[0]= mvConNode[0]->getX() - mvConNode[1]->getX();
            vP1P0[1]= mvConNode[0]->getY() - mvConNode[1]->getY();
            vP1P0[2]= mvConNode[0]->getZ() - mvConNode[1]->getZ();

            //正規化
            Normalized(vPOP);  Normalized(vP2P3);  Normalized(vP1P0);

            //-----------------------------------------------------//
            //vPOPとp2-p3ベクトルの最近接点 (以下の３つのルーチンはセットで使用)
            //
            NearCrossP(vPOP,vP2P3, vOP,mvConNode[2]);
            s_NearCrossP();
            t_NearCrossP();

            //vPOPベクトル,P2P3ベクトルの最近接点2点 NP0,NP1
            vNP0[0]= mS*vPOP[0] + vOP[0];//ベクトルvPOP * s + ベクトル基点
            vNP0[1]= mS*vPOP[1] + vOP[1];
            vNP0[2]= mS*vPOP[2] + vOP[2];

            vNP1[0]= mT*vP2P3[0] + mvConNode[2]->getX();//ベクトルvP2P3 * t + ベクトル基点
            vNP1[1]= mT*vP2P3[1] + mvConNode[2]->getY();
            vNP1[2]= mT*vP2P3[2] + mvConNode[2]->getZ();

            //vPOPとvP2P3の近似交点
            vInterP0[0]= (vNP0[0]+vNP1[0])*0.5;
            vInterP0[1]= (vNP0[1]+vNP1[1])*0.5;
            vInterP0[2]= (vNP0[2]+vNP1[2])*0.5;
                
            //InterP0の補間値
            //
            dLenP2P3= VectorLength(mvConNode[2],mvConNode[3]);
            dLenP2InP= VectorLength(mvConNode[2],vInterP0);
            //Valueをそれぞれ線形補間
            for(idof=0; idof < numOfDOF; idof++){
                vValInterP0[idof]= LinearInter(dLenP2InP/dLenP2P3,idof,mvConNode[2],mvConNode[3],valType);
            };
            //-----------------------------------------------------//


            //-----------------------------------------------------//
            //vPOPとp1-P0ベクトルの最近接点 (以下の３つのルーチンはセットで使用)
            //--
            NearCrossP(vPOP,vP1P0, vOP,mvConNode[1]);
            s_NearCrossP();
            t_NearCrossP();

            //vPOPベクトル,P2P3ベクトルの最近接点2点 NP0,NP1
            vNP0[0]= mS*vPOP[0] + vOP[0];//ベクトルvPOP * s + ベクトル基点
            vNP0[1]= mS*vPOP[1] + vOP[1];
            vNP0[2]= mS*vPOP[2] + vOP[2];

            vNP1[0]= mT*vP1P0[0] + mvConNode[1]->getX();//ベクトルvP1P0* t + ベクトル基点
            vNP1[1]= mT*vP1P0[1] + mvConNode[1]->getY();
            vNP1[2]= mT*vP1P0[2] + mvConNode[1]->getZ();

            //vPOPとvP1P0の近似交点
            vInterP1[0]= (vNP0[0]+vNP1[0])*0.5;
            vInterP1[1]= (vNP0[1]+vNP1[1])*0.5;
            vInterP1[2]= (vNP0[2]+vNP1[2])*0.5;

            //InterP1の補間値
            //
            dLenP0P1= VectorLength(mvConNode[1],mvConNode[0]);
            dLenP1InP= VectorLength(mvConNode[1],vInterP1);
            //Valueをそれぞれ線形補間
            for(idof=0; idof < numOfDOF; idof++){
                vValInterP1[idof]= LinearInter(dLenP1InP/dLenP0P1,idof,mvConNode[1],mvConNode[0],valType);
            };
            //-----------------------------------------------------//
            

            // Coefの計算
            // --
            // Quad頂点0 のCoef
            mvvCoef[islave][0]= CoefTermA(mvConNode[1],vInterP1,mvConNode[0],mvConNode[1]) * CoefTermB(vInterP0,pSlaveP, vInterP1,vInterP0);
            // Quad頂点1 のCoef
            mvvCoef[islave][1]= CoefTermA(mvConNode[0],vInterP1,mvConNode[0],mvConNode[1]) * CoefTermB(vInterP0,pSlaveP, vInterP1,vInterP0);
            // Quad頂点2 のCoef
            mvvCoef[islave][2]= CoefTermA(mvConNode[3],vInterP0,mvConNode[2],mvConNode[3]) * CoefTermB(vInterP1,pSlaveP, vInterP1,vInterP0);
            // Quad頂点3 のCoef
            mvvCoef[islave][3]= CoefTermA(mvConNode[2],vInterP0,mvConNode[2],mvConNode[3]) * CoefTermB(vInterP1,pSlaveP, vInterP1,vInterP0);
            
            // 2次要素の場合
            if(mnOrder==ElementOrder::Second){
                //    mvvCoef[islave][0] *= 0.5;
                //    mvvCoef[islave][1] *= 0.5;
                //    mvvCoef[islave][2] *= 0.5;
                //    mvvCoef[islave][3] *= 0.5;
                //    mvvCoef[islave][4] = (mvvCoef[islave][0] + mvvCoef[islave][1])*0.5;
                //    mvvCoef[islave][5] = (mvvCoef[islave][1] + mvvCoef[islave][2])*0.5;
                //    mvvCoef[islave][6] = (mvvCoef[islave][2] + mvvCoef[islave][3])*0.5;
                //    mvvCoef[islave][7] = (mvvCoef[islave][3] + mvvCoef[islave][0])*0.5;
                
                // 頂点のCoefだけ有効 : 辺の節点Ceof=0.0  2011.01.13
                //
                mvvCoef[islave][4] = 0.0;
                mvvCoef[islave][5] = 0.0;
                mvvCoef[islave][6] = 0.0;
                mvvCoef[islave][7] = 0.0;
            }
            
            break;

        case(3)://三角形
            //E-F ベクトルの交点==p2
            vOP[0]= mvConNode[2]->getX();
            vOP[1]= mvConNode[2]->getY();
            vOP[2]= mvConNode[2]->getZ();

            //OPとスレーブ点のベクトル
            vPOP[0]= pSlaveP->getX() - vOP[0];
            vPOP[1]= pSlaveP->getY() - vOP[1];
            vPOP[2]= pSlaveP->getZ() - vOP[2];
            //p1-p0ベクトル
            vP1P0[0]= mvConNode[0]->getX() - mvConNode[1]->getX();
            vP1P0[1]= mvConNode[0]->getY() - mvConNode[1]->getY();
            vP1P0[2]= mvConNode[0]->getZ() - mvConNode[1]->getZ();

            //正規化
            Normalized(vPOP);  Normalized(vP1P0);

            //vPOPとp1-P0ベクトルの最近接点 (以下の３つのルーチンはセットで使用)
            //--
            NearCrossP(vPOP,vP1P0, vOP,mvConNode[1]);
            s_NearCrossP();
            t_NearCrossP();

            //vPOPベクトル,P2P3ベクトルの最近接点2点 NP0,NP1
            vNP0[0]= mS*vPOP[0] + vOP[0];//ベクトルvPOP * s + ベクトル基点
            vNP0[1]= mS*vPOP[1] + vOP[1];
            vNP0[2]= mS*vPOP[2] + vOP[2];

            vNP1[0]= mT*vP1P0[0] + mvConNode[1]->getX();//ベクトルvP1P0* t + ベクトル基点
            vNP1[1]= mT*vP1P0[1] + mvConNode[1]->getY();
            vNP1[2]= mT*vP1P0[2] + mvConNode[1]->getZ();

            //vPOPとvP1P0の近似交点
            vInterP1[0]= (vNP0[0]+vNP1[0])*0.5;
            vInterP1[1]= (vNP0[1]+vNP1[1])*0.5;
            vInterP1[2]= (vNP0[2]+vNP1[2])*0.5;
            //InterP1の補間値
            //
            dLenP0P1= VectorLength(mvConNode[1],mvConNode[0]);
            dLenP1InP= VectorLength(mvConNode[1],vInterP1);
            //Valueを線形補間
            for(idof=0; idof < numOfDOF; idof++){
                vValInterP1[idof]= LinearInter(dLenP1InP/dLenP0P1,idof,mvConNode[1],mvConNode[0],valType);
            };


            //vInterP0 == p2
            vInterP0[0]= mvConNode[2]->getX();
            vInterP0[1]= mvConNode[2]->getY();
            vInterP0[2]= mvConNode[2]->getZ();
            //vValInterP0
            switch(valType){
                case(MPCValueType::Displacement):
                    //x,y,z の変位値
                    for(idof=0; idof < numOfDOF; idof++) vValInterP0[idof]= mvConNode[2]->getDisp(idof);
                    break;
                case(MPCValueType::Scalar):
                    //スカラー値
                    for(idof=0; idof < numOfDOF; idof++) vValInterP0[idof]= mvConNode[2]->getScalar(idof);
                    break;
            }
            
            // Coefの計算 :: 三角形の面積比に改訂
            //   # CrossProduct()は,mvVectorに計算結果が入る.
            // 
            CrossProduct(mvConNode[0], 1,2);   //三角形全体の外積
            aArea= VectorLength(mvVector)*0.5;//mvVectorの大きさの1/2=三角形面積
            CrossProduct(pSlaveP, 1,2);
            vArea[0]= VectorLength(mvVector)*0.5;
            CrossProduct(pSlaveP, 2,0);
            vArea[1]= VectorLength(mvVector)*0.5;
            CrossProduct(pSlaveP, 0,1);
            vArea[2]= VectorLength(mvVector)*0.5;

            //// slaveID単位のCoef準備
            //vuint vVertCoef; vVertCoef.resize(3);
            //mmvCoef[slaveID]= vVertCoef;

            mvvCoef[islave][0]= vArea[0]/aArea;// Tri 頂点0 のCoef
            mvvCoef[islave][1]= vArea[1]/aArea;// Tri 頂点1 のCoef
            mvvCoef[islave][2]= vArea[2]/aArea;// Tri 頂点2 のCoef

            // 2次要素の場合
            if(mnOrder==ElementOrder::Second){
                //    mvvCoef[islave][0] *= 0.5;
                //    mvvCoef[islave][1] *= 0.5;
                //    mvvCoef[islave][2] *= 0.5;
                //    mvvCoef[islave][3] = (mvvCoef[islave][0] + mvvCoef[islave][1])*0.5;
                //    mvvCoef[islave][4] = (mvvCoef[islave][1] + mvvCoef[islave][2])*0.5;
                //    mvvCoef[islave][5] = (mvvCoef[islave][2] + mvvCoef[islave][0])*0.5;

                // 頂点のCoefだけ有効 : 辺の節点Ceof=0.0  2011.01.13
                //
                mvvCoef[islave][3] = 0.0;
                mvvCoef[islave][4] = 0.0;
                mvvCoef[islave][5] = 0.0;
            }

            break;

        case(2)://Beam(線要素)
            vInterP0[0]= mvConNode[0]->getX();
            vInterP0[1]= mvConNode[0]->getY();
            vInterP0[2]= mvConNode[0]->getZ();

            vInterP1[0]= mvConNode[1]->getX();
            vInterP1[1]= mvConNode[1]->getY();
            vInterP1[2]= mvConNode[1]->getZ();

            switch(valType){
                case(MPCValueType::Displacement):
                    for(idof=0; idof < numOfDOF; idof++){
                        vValInterP0[idof]= mvConNode[0]->getDisp(idof);
                        vValInterP1[idof]= mvConNode[1]->getDisp(idof);
                    };
                    break;
                case(MPCValueType::Scalar):
                    for(idof=0; idof < numOfDOF; idof++){
                        vValInterP0[idof]= mvConNode[0]->getScalar(idof);
                        vValInterP1[idof]= mvConNode[1]->getScalar(idof);
                    };
                    break;
            }
            //// slaveID単位のCoef準備
            //vuint vVertCoef; vVertCoef.resize(2);
            //mmvCoef[slaveID]= vVertCoef;

            // Coefの計算
            //
            mvvCoef[islave][0]= CoefTermB(vInterP1,pSlaveP, vInterP0,vInterP1);// Beam 頂点0 のCoef
            mvvCoef[islave][1]= CoefTermB(vInterP0,pSlaveP, vInterP0,vInterP1);// Beam 頂点1 のCoef

            // 2次要素の場合
            if(mnOrder==ElementOrder::Second){
                //mvvCoef[islave][0] *= 0.5;
                //mvvCoef[islave][1] *= 0.5;
                //mvvCoef[islave][2] = (mvvCoef[islave][0] + mvvCoef[islave][1])*0.5;

                // 頂点のCoefだけ有効 : 辺の節点Ceof=0.0  2011.01.13
                //
                mvvCoef[islave][2] = 0.0;
            }

            break;

        default:
            break;
    }
    // pSlaveP:                 スレーブ点(このルーチンの最初に取得ずみ.)
    // vInterP0,vInterP1:       vPOPベクトルと辺の交点座標
    // vValInterP0,vValInterP1: 辺上のDispの補間値
    // --
    // スレーブ点のDispを補間
    // --
    dLenInP0InP1= VectorLength(vInterP0, vInterP1); //InterP0-InterP1の距離
    dLenInP0SlaveP= VectorLength(vInterP0, pSlaveP);//InterP0-SlavePの距離
    // Valueをそれぞれ線形補間
    for(idof=0; idof < numOfDOF; idof++){
        vValSlaveP[idof]= LinearInter(dLenInP0SlaveP/dLenInP0InP1,idof, vValInterP0,vValInterP1);//スレーブ点での補間値
    };

    // 補間された値をスレーブ点に代入
    //
    switch(valType){
        case(MPCValueType::Displacement):
            for(idof=0; idof < numOfDOF; idof++)
                pSlaveP->setDisp(idof, vValSlaveP[idof]);
            break;
        case(MPCValueType::Scalar):
            for(idof=0; idof < numOfDOF; idof++)
                pSlaveP->setScalar(idof, vValSlaveP[idof]);
            break;
    }
}


// スレーブ点ID,頂点番号を引数としてCoef を提供  => ConNodeからCoefを取得させる方向に変更.'10.02.23
//   Coef== mvvCoef[islave][ivert]
double& CMasterFace::getCoef(const uiint& slaveID, const uiint& ivert)
{
    uiint islave= mmSlaveID2Index[slaveID];

    return mvvCoef[islave][ivert];
}













