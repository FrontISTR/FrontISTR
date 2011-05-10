//
//  CommMesh2.cpp
//
//
//
//                  2010.03.02
//                  k.Takeda
#include <vector>

#include "CommMesh2.h"
using namespace pmw;

CCommMesh2::CCommMesh2()
{
    ;
}
CCommMesh2::~CCommMesh2()
{
    uint numOfNode= mvCommNode.size();
    CCommNode *pCommNode;
    uint inode;
    for(inode=0; inode< numOfNode; inode++){
        pCommNode= mvCommNode[inode];

        if(pCommNode->getLevel()==mMGLevel){
            delete pCommNode;
        }
    };

    uint numOfFace= mvCommFace.size();
    CCommFace *pCommFace;
    uint iface;
    for(iface=0; iface< numOfFace; iface++){
        pCommFace= mvCommFace[iface];
        delete pCommFace;
    };
}

void CCommMesh2::addCommNode(CCommNode* pCommNode)
{
    mvCommNode.push_back(pCommNode);
    mmCommNodeID2Index[pCommNode->getID()]= mvCommNode.size()-1;
}
void CCommMesh2::addCommFace(CCommFace* pCommFace)
{
    mvCommFace.push_back(pCommFace);
    mmCommFaceID2Index[pCommFace->getID()]= mvCommFace.size()-1;
}
CCommNode* CCommMesh2::getCommNode(const uint& id)
{
    uint index;
    index= mmCommNodeID2Index[id];

    return mvCommNode[index];
}
CCommFace* CCommMesh2::getCommFace(const uint& id)
{
    uint index;
    index= mmCommFaceID2Index[id];

    return mvCommFace[index];
}


////// 頂点CommNodeのProgCommMesh2への代入
//////
////void CCommMesh2::setupVertCommNode(CCommMesh2* pProgCommMesh)
//
// ---- 
// ProgCommMesh2 へ mvCommNodeを代入
// ----
void CCommMesh2::setupCommNode(CCommMesh2* pProgCommMesh)
{
    uint numOfNode;
    CCommNode *pCommNode;
    
    numOfNode= mvCommNode.size();
    pProgCommMesh->reserveCommNode(numOfNode);

    uint inode;
    for(inode=0; inode< numOfNode; inode++){
        pCommNode= mvCommNode[inode];
        pProgCommMesh->addCommNode(pCommNode);
    };
}


// 頂点へのCommFace集合
//
void CCommMesh2::setupAggFace()
{
    CCommFace *pCommFace;
    CCommNode *pCommNode;

    uint numOfCommFace(mvCommFace.size());
    uint numOfCommNode;

    uint iface,inode;

    //------------------------------
    // 第2段目以降のAggFaceのために
    //   下位Levelでセットされた,FaceIDをクリア
    //------------------------------
    numOfCommNode= mvCommNode.size();
    for(inode=0; inode< numOfCommNode; inode++){
        pCommNode= mvCommNode[inode];

        pCommNode->clearAggElemID();
        pCommNode->clearNeibElemVert();
    };

    

    for(iface=0; iface< numOfCommFace; iface++){
        pCommFace= mvCommFace[iface];

        uint nNumOfVert= pCommFace->getNumOfVert();//辺ノードがセットされる前なので頂点のみ
        
        for(inode=0; inode< nNumOfVert; inode++){
            pCommNode= pCommFace->getCommNode(inode);

            pCommNode->setAggElemID(pCommFace->getID());
            pCommNode->setNeibElemVert(pCommFace->getID(), inode);
        };
    };
    
}
// 辺への新CommNodeの追加
//
void CCommMesh2::setupEdgeCommNode(CCommMesh2 *pProgCommMesh, const uint& nLevel)
{
    PairCommNode pairCommNode;
    CCommFace *pFace, *pNeibFace;
    CCommNode* pEdgeCommNode;
    vdouble vCoord; vCoord.resize(3);
    uint numOfFace, numOfEdge;
    uint iface,iedge;

    uint countID= mvCommNode.size();//辺ノードID初期値
    mEdgeCount= 0;//辺に生成するノード数のカウンター

    numOfFace= mvCommFace.size();
    for(iface=0; iface< numOfFace; iface++){
        
        pFace= mvCommFace[iface];
        
        // 辺ノード生成 実行判定
        bool bExec(false);
        if(nLevel > 0 && pFace->getOrder()==ElementOrder::Second) bExec=true;
        if(pFace->getOrder()==ElementOrder::First) bExec=true;


        // 辺ノード生成
        //
        if(bExec){
            numOfEdge= pFace->getNumOfEdge();

            for(iedge=0; iedge< numOfEdge; iedge++){
                if(!pFace->isEdgeNodeMarking(iedge)){
                    pairCommNode= pFace->getEdgePairCommNode(iedge);

                    uint numOfAggFaceA= pairCommNode.first->getNumOfAggElem();
                    uint numOfAggFaceB= pairCommNode.second->getNumOfAggElem();
                    uint faceIDa,faceIDb;
                    uint iagg,jagg;
                    bool bfind(false);//隣接Faceが見つかった場合
                    // ----
                    // 隣接するFaceを検索してセット .&&. 辺CommNodeを生成してセット
                    // ----
                    for(iagg=0; iagg< numOfAggFaceA; iagg++){
                        faceIDa= pairCommNode.first->getAggElemID(iagg);

                        for(jagg=0; jagg< numOfAggFaceB; jagg++){
                            faceIDb= pairCommNode.second->getAggElemID(jagg);

                            // 辺に接続するFace(pFace自身は除外)
                            if(faceIDa==faceIDb && faceIDa!=pFace->getID()){

                                // 隣のCommFace "ID => Index" 2010.05.12
                                uint index= mmCommFaceID2Index[faceIDa];
                                pNeibFace= mvCommFace[index];

                                pEdgeCommNode= new CCommNode;// <<<<- CommNode生成

                                //------------------------
                                // 辺ノードへの座標,IDのセット
                                //------------------------
                                vCoord[0]= pairCommNode.first->getX(); vCoord[0]+= pairCommNode.second->getX(); vCoord[0] *= 0.5;
                                vCoord[1]= pairCommNode.first->getY(); vCoord[1]+= pairCommNode.second->getY(); vCoord[1] *= 0.5;
                                vCoord[2]= pairCommNode.first->getZ(); vCoord[2]+= pairCommNode.second->getZ(); vCoord[2] *= 0.5;
                                pEdgeCommNode->setCoord(vCoord);

                                pEdgeCommNode->setLevel(mMGLevel+1);//2010.05.27
                                pEdgeCommNode->setID(countID);
                                countID++;//////////  ID ,次の為にカウントアップ
                                mEdgeCount++;// FaceNodeの為にカウントアップ

                                if(pProgCommMesh)
                                    pProgCommMesh->addCommNode(pEdgeCommNode);// 上位CommMesh2 にノードを追加

                                // 2次要素面の場合
                                //
                                if(pFace->getOrder()==ElementOrder::Second){
                                    mvCommNode.push_back(pEdgeCommNode);// 自身のノード集合に追加
                                }


                                // pFaceへの処理
                                pFace->setEdgeCommFace(pNeibFace, iedge);
                                pFace->setEdgeCommNode(pEdgeCommNode, iedge);
                                pFace->markingEdgeNode(iedge);


                                // pNeibFaceへの処理
                                pNeibFace->setEdgeCommFace(pFace, pairCommNode);
                                pNeibFace->setEdgeCommNode(pEdgeCommNode, pairCommNode);
                                pNeibFace->markingEdgeNode(pairCommNode);


                                bfind=true;//1辺に隣接するFaceは一つなので(pFace自身は除外してる)
                                break;
                            }
                        };// jaggループ

                        if(bfind) break;//iaggのループをブレイク(既に隣接Faceを発見したので)

                    };// iaggループ


                    // 隣接していない辺にCommNodeを生成
                    // ----
                    if(!bfind){

                        pEdgeCommNode= new CCommNode;// <<<<<<<<<- CommNode生成

                        //------------------------
                        // 辺ノードへの座標,IDのセット
                        //------------------------
                        vCoord[0]= pairCommNode.first->getX(); vCoord[0]+= pairCommNode.second->getX(); vCoord[0] *= 0.5;
                        vCoord[1]= pairCommNode.first->getY(); vCoord[1]+= pairCommNode.second->getY(); vCoord[1] *= 0.5;
                        vCoord[2]= pairCommNode.first->getZ(); vCoord[2]+= pairCommNode.second->getZ(); vCoord[2] *= 0.5;
                        pEdgeCommNode->setCoord(vCoord);

                        pEdgeCommNode->setLevel(mMGLevel+1);//2010.05.27
                        pEdgeCommNode->setID(countID);
                        countID++;////////// ID ,次の為にカウントアップ
                        mEdgeCount++;// FaceNodeの為にカウントアップ

                        if(pProgCommMesh)
                            pProgCommMesh->addCommNode(pEdgeCommNode);// 上位CommMesh2にノードを追加

                        // 2次要素面の場合
                        //
                        if(pFace->getOrder()==ElementOrder::Second){
                            mvCommNode.push_back(pEdgeCommNode);// 自身のノード集合に追加
                        }

                        // pFaceへの処理
                        pFace->setEdgeCommNode(pEdgeCommNode, iedge);
                        pFace->markingEdgeNode(iedge);
                    }
                }// !isEdgeNodeMarking
            };// iedgeループ
            
            pFace->replaceEdgeCommNode();// 2次ノード対応
        
        }// if(bExec)
        
    };// ifaceループ

}

// 面中心への新CommNodeの追加
//
void CCommMesh2::setupFaceCommNode(CCommMesh2 *pProgCommMesh)
{
    uint level;
    CCommNode *pFaceCommNode;
    
    CCommFace *pFace;
    uint numOfFace= mvCommFace.size();
    uint iface;
    //uint countID= mvCommNode.size()+mvEdgeCommNode.size();
    uint countID= mvCommNode.size() + mEdgeCount;

    for(iface=0; iface< numOfFace; iface++){
        pFace= mvCommFace[iface];
        
        //-- Quad,TriangleのみにFaceNodeを追加 --//
        if(pFace->getNumOfEdge() > 2){
            pFaceCommNode= new CCommNode;////////// <<<--- FaceのCommNode

            pFace->setFaceCommNode(pFaceCommNode);

            vdouble vCoord; vCoord.resize(3);
            uint numOfVert= pFace->getNumOfVert();
            CCommNode *pVertCommNode;
            vCoord[0]=0.0; vCoord[1]=0.0; vCoord[2]=0.0;
            uint ivert;
            level=0;
            for(ivert=0; ivert< numOfVert; ivert++){
                pVertCommNode= pFace->getCommNode(ivert);

                ////CommNodeの中から最大のLevelを見つける
                //if(pVertCommNode->getLevel() > level) level= pVertCommNode->getLevel();

                vCoord[0] += pVertCommNode->getX();  vCoord[1] += pVertCommNode->getY();  vCoord[2] += pVertCommNode->getZ();
            };
            vCoord[0] /= numOfVert; vCoord[1] /= numOfVert; vCoord[2] /= numOfVert;

            //level++;//面を構成しているCommNodeより一段上へ.

            pFaceCommNode->setLevel(mMGLevel+1);//2010.05.27
            pFaceCommNode->setCoord(vCoord);
            pFaceCommNode->setID(countID);
            countID++;////////// ID カウントアップ

            pProgCommMesh->addCommNode(pFaceCommNode);
            //// mvFaceCommNode.push_back(pFaceCommNode);
        }
    };
}


// Refine後のvector解放(辺ノード)
//
void CCommMesh2::deleteProgData()
{
    uint iface, nNumOfFace = mvCommFace.size();

    for(iface=0; iface < nNumOfFace; iface++) mvCommFace[iface]->deleteProgData();
}









