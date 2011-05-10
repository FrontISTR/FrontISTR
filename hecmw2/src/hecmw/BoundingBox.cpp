//
//  BoundingBox.cpp
//
//
//
//                  2009.10.14
//                  2009.10.14
//                  k.Takeda
#include "Vertex.h"
#include "BoundingBox.h"
using namespace pmw;

// construct & destruct
//
CBoundingBox::CBoundingBox()
{
    mMinCoord.resize(3);
    mMaxCoord.resize(3);

    mvOBBCenter.resize(3);//OBB座標の中心点(x,y,z)
    mvOBBX.resize(3);//OBB座標のX方向ベクトル
    mvOBBY.resize(3);//OBB座標のY方向ベクトル
    mvOBBZ.resize(3);//OBB座標のZ方向ベクトル
    mvE.resize(3);//OBB座標のX,Y,Z方向の範囲
}

CBoundingBox::~CBoundingBox()
{
}

//void CBoundingBox::sizing()
//{
//    mvPoint[0]->setX(mMinCoord[0]); mvPoint[0]->setY(mMinCoord[1]); mvPoint[0]->setZ(mMinCoord[2]);//mPoint[0]:最小
//    mvPoint[1]->setX(mMaxCoord[0]); mvPoint[1]->setY(mMinCoord[1]); mvPoint[1]->setZ(mMinCoord[2]);
//    mvPoint[2]->setX(mMaxCoord[0]); mvPoint[2]->setY(mMinCoord[1]); mvPoint[2]->setZ(mMaxCoord[2]);
//    mvPoint[3]->setX(mMinCoord[0]); mvPoint[3]->setY(mMinCoord[1]); mvPoint[3]->setZ(mMaxCoord[2]);
//    mvPoint[4]->setX(mMinCoord[0]); mvPoint[4]->setY(mMaxCoord[1]); mvPoint[4]->setZ(mMinCoord[2]);
//    mvPoint[5]->setX(mMaxCoord[0]); mvPoint[5]->setY(mMaxCoord[1]); mvPoint[5]->setZ(mMinCoord[2]);
//    mvPoint[6]->setX(mMaxCoord[0]); mvPoint[6]->setY(mMaxCoord[1]); mvPoint[6]->setZ(mMaxCoord[2]);//mPoint[6]:最大
//    mvPoint[7]->setX(mMinCoord[0]); mvPoint[7]->setY(mMaxCoord[1]); mvPoint[7]->setZ(mMaxCoord[2]);
//}

// ABBのサイズ決定
// 
//  BBを作るマスター面の最大,最小を取得して,BBのサイズ(座標の最大,最小)を決める.
//
void CBoundingBox::sizingABB(CSkinFace* pFace)
{
    uint numOfVert= pFace->getNumOfNode();
    uint ivert;
    CContactNode* pConNode;
    // 初期値
    pConNode= pFace->getNode(0);
    mMinCoord[0]= pConNode->getX(); mMinCoord[1]= pConNode->getY(); mMinCoord[2]= pConNode->getZ();
    mMaxCoord[0]= pConNode->getX(); mMaxCoord[1]= pConNode->getY(); mMaxCoord[2]= pConNode->getZ();
    // 初期値以外と比較
    for(ivert=1; ivert< numOfVert; ivert++){
        pConNode= pFace->getNode(ivert);

        if(mMinCoord[0] > pConNode->getX()) mMinCoord[0]= pConNode->getX();
        if(mMinCoord[1] > pConNode->getY()) mMinCoord[1]= pConNode->getY();
        if(mMinCoord[2] > pConNode->getZ()) mMinCoord[2]= pConNode->getZ();

        if(mMaxCoord[0] < pConNode->getX()) mMaxCoord[0]= pConNode->getX();
        if(mMaxCoord[1] < pConNode->getY()) mMaxCoord[1]= pConNode->getY();
        if(mMaxCoord[2] < pConNode->getZ()) mMaxCoord[2]= pConNode->getZ();

    };//ivertループ
}


// --
// ABB判定 Face(スレーブ面)
// --
bool CBoundingBox::judgeABB(CSkinFace* pFace)
{
    CContactNode* pConNode;
    uint numOfConNode= pFace->getNumOfNode();
    
    uint icnode;
    bool bCheck(true);//初期値：TRUE
    for(icnode=0; icnode< numOfConNode; icnode++){
        pConNode= pFace->getNode(icnode);

        //一つでもFALSEならば,FaceはFALSE
        if(!judgeABB(pConNode)){
            bCheck=false;
            break;
        }
    };
    
    return bCheck;
}
// --
// ABB判定 ContactNode(スレーブConNode)
// --
bool CBoundingBox::judgeABB(CContactNode* pConNode)
{
    bool bCheck(false);

    double X,Y,Z;
    X=pConNode->getX(); Y=pConNode->getY(); Z=pConNode->getZ();

    if(X >= mMinCoord[0] && X <= mMaxCoord[0]){
        if(Y >= mMinCoord[1] && Y <= mMaxCoord[1]){
            if(Z >= mMinCoord[2] && Z <= mMaxCoord[2]){
                bCheck=true;
            }
        }
    }

    return bCheck;
}



// OBB座標軸 && 範囲 の設定
//
//  マスター面は,必ず"面"であることを前提としている. -> Beamによるマスター面？は,まだ考慮していない.
//  -----------------------------------------
//
void CBoundingBox::sizingOBB(CSkinFace* pFace)
{
    //OBB 中心点
    uint numOfVert= pFace->getNumOfNode();
    uint ivert;
    CContactNode* pConNode;
    mvOBBCenter[0]=0.0; mvOBBCenter[1]=0.0; mvOBBCenter[2]=0.0;
    for(ivert=0; ivert< numOfVert; ivert++){
        pConNode= pFace->getNode(ivert);

        mvOBBCenter[0]+= pConNode->getX();
        mvOBBCenter[1]+= pConNode->getY();
        mvOBBCenter[2]+= pConNode->getZ();
    };
    mvOBBCenter[0]/= (double)numOfVert;
    mvOBBCenter[1]/= (double)numOfVert;
    mvOBBCenter[2]/= (double)numOfVert;

    //    //debug
    //    cout << "BoundingBox::center_X= " << mvOBBCenter[0]
    //                    << ", center_y= " << mvOBBCenter[1]
    //                    << ", center_z= " << mvOBBCenter[2] << endl;


    //OBB座標 Z方向 (3つの頂点による外積) :: 面であることが前提の処理(Beamは,除外)
    //
    CContactNode *pVert0,*pVert1,*pVert2;
    pVert0= pFace->getNode(0); pVert1= pFace->getNode(1); pVert2= pFace->getNode(2);
    
    // 外積に使用する2つのベクトル
    double X1= pVert1->getX() - pVert0->getX();
    double Y1= pVert1->getY() - pVert0->getY();
    double Z1= pVert1->getZ() - pVert0->getZ();

    double X2= pVert2->getX() - pVert0->getX();
    double Y2= pVert2->getY() - pVert0->getY();
    double Z2= pVert2->getZ() - pVert0->getZ();

    // 外積
    // x = y1 z2 - z1 y2
    // y = z1 x2 - x1 z2
    // z = x1 y2 - y1 x2
    //
    mvOBBZ[0]= Y1*Z2 - Z1*Y2;
    mvOBBZ[1]= Z1*X2 - X1*Z2;
    mvOBBZ[2]= X1*Y2 - Y1*X2;

    // 正規化
    double absValue= sqrt(mvOBBZ[0]*mvOBBZ[0] + mvOBBZ[1]*mvOBBZ[1] + mvOBBZ[2]*mvOBBZ[2]);

    mvOBBZ[0] /= absValue;
    mvOBBZ[1] /= absValue;
    mvOBBZ[2] /= absValue;

    //    //debug
    //    cout << "BoundingBox::obbZ_X= " << mvOBBZ[0]
    //                    << ", obbZ_y= " << mvOBBZ[1]
    //                    << ", obbZ_z= " << mvOBBZ[2] << endl;
    

    //OBB座標 X方向(頂点0と頂点1の中間点)
    mvOBBX[0]= (pVert0->getX() + pVert1->getX())*0.5;
    mvOBBX[1]= (pVert0->getY() + pVert1->getY())*0.5;
    mvOBBX[2]= (pVert0->getZ() + pVert1->getZ())*0.5;
    //ベクトル化
    mvOBBX[0] -= mvOBBCenter[0];
    mvOBBX[1] -= mvOBBCenter[1];
    mvOBBX[2] -= mvOBBCenter[2];

    // 正規化
    absValue= sqrt(mvOBBX[0]*mvOBBX[0] + mvOBBX[1]*mvOBBX[1] + mvOBBX[2]*mvOBBX[2]);

    mvOBBX[0] /= absValue;
    mvOBBX[1] /= absValue;
    mvOBBX[2] /= absValue;

    //    //debug
    //    cout << "BoundingBox::obbX_X= " << mvOBBX[0]
    //                    << ", obbX_y= " << mvOBBX[1]
    //                    << ", obbX_z= " << mvOBBX[2] << endl;

    //OBB座標 Y方向
    X1= mvOBBZ[0];
    Y1= mvOBBZ[1];
    Z1= mvOBBZ[2];

    X2= mvOBBX[0];
    Y2= mvOBBX[1];
    Z2= mvOBBX[2];

    mvOBBY[0]= Y1*Z2 - Z1*Y2;
    mvOBBY[1]= Z1*X2 - X1*Z2;
    mvOBBY[2]= X1*Y2 - Y1*X2;

    // 正規化
    absValue= sqrt(mvOBBY[0]*mvOBBY[0] + mvOBBY[1]*mvOBBY[1] + mvOBBY[2]*mvOBBY[2]);

    mvOBBY[0] /= absValue;
    mvOBBY[1] /= absValue;
    mvOBBY[2] /= absValue;

    //    //debug
    //    cout << "BoundingBox::obbY_X= " << mvOBBY[0]
    //                    << ", obbY_y= " << mvOBBY[1]
    //                    << ", obbY_z= " << mvOBBY[2] << endl;


    //範囲(Ex,Ey,Ez) :: mvE[0],mvE[1],mvE[2]
    //
    mvE[0]=0.0; mvE[1]=0.0; mvE[2]=0.0;//OBBX,OBBY,OBBZの範囲値クリア
    double obbX,obbY;//tempo変数
    
    for(ivert=0; ivert< numOfVert; ivert++){
        pConNode= pFace->getNode(ivert);

        X1= pConNode->getX() - mvOBBCenter[0];
        Y1= pConNode->getY() - mvOBBCenter[1];
        Z1= pConNode->getZ() - mvOBBCenter[2];

        // OBB_X,OBB_Y 座標軸への投影値から最大値を取得
        // --
        obbX= abs(X1*mvOBBX[0] + Y1*mvOBBX[1] + Z1*mvOBBX[2]);//OBB_X への投影(内積)の絶対値
        if(ivert==0) mvE[0]=obbX;//初期値
        if(obbX > mvE[0]) mvE[0]=obbX;//OBB_X への投影値の絶対値が最大になるサイズをExとする.
        
        obbY= abs(X1*mvOBBY[0] + Y1*mvOBBY[1] + Z1*mvOBBY[2]);//OBB_Y への投影(内積)の絶対値
        if(ivert==0) mvE[1]=obbY;//初期値
        if(obbY > mvE[1]) mvE[1]=obbY;//OBB_Y への投影値の絶対値が最大になるサイズをEyとする.


        // OBB_Z 座標軸 (面ベクトル方向)
        // --
        mvE[2]=(mvE[0]+mvE[1])*0.01;//面ベクトルなので,厚さ==0.0 -> 閾値として(Ex+Ey)*0.01 を使用
    };

    //範囲(Ex,Ey,Ez)を拡大(10%)
    mvE[0] *= 1.10;
    mvE[1] *= 1.10;
    mvE[2] *= 1.10;

    //    //debug
    //    cout << "BoundingBox::Ex= " << mvE[0] << ", Ey= " << mvE[1] << ", Ez= " << mvE[2] << endl;
}

// --
// OBB判定 ContactNode(スレーブConNode)
// --
bool CBoundingBox::judgeOBB(CContactNode* pConNode)
{
    bool bCheck(false);
    
    // (スレーブ点-OBB中心)ベクトル
    double X1= pConNode->getX() - mvOBBCenter[0];
    double Y1= pConNode->getY() - mvOBBCenter[1];
    double Z1= pConNode->getZ() - mvOBBCenter[2];

    double projX,projY,projZ;//projection OBB coord
    // 各OBB座標軸への投影(内積):スカラー値
    // 
    projX= abs(X1*mvOBBX[0] + Y1*mvOBBX[1] + Z1*mvOBBX[2]);//OBB_X への投影の絶対値
    projY= abs(X1*mvOBBY[0] + Y1*mvOBBY[1] + Z1*mvOBBY[2]);//OBB_Y への投影の絶対値
    projZ= abs(X1*mvOBBZ[0] + Y1*mvOBBZ[1] + Z1*mvOBBZ[2]);//OBB_Z への投影の絶対値

    
    // Ex,Ey,Ezとの比較(BoundingBoxに入っているか?)
    //  *絶対値による比較*
    if(projX <= mvE[0]){
        if(projY <= mvE[1]){
            if(projZ <= mvE[2]) bCheck=true;
        }
    }
    
    //    //debug
    //    if(bCheck && pConNode->getID()==22){
    //        cout << "conID==22" << endl;
    //        cout << "BoundingBox::Ex= " << mvE[0]
    //                        << ", Ey= " << mvE[1]
    //                        << ", Ez= " << mvE[2]  << endl;
    //        cout << "BoundingBox::judgeOBB, projX= " << projX
    //                                  << ", projY= " << projY
    //                                  << ", projZ= " << projZ << endl;
    //    }

    return bCheck;
}
































