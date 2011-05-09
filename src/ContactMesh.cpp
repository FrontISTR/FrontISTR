
#include <vector>
#include "Vertex.h"
#include "ContactNode.h"

//
// ContactMesh.cpp
//
//			2009.01.08
//			2009.01.08
//			k.Takeda
#include "ContactMesh.h"
using namespace pmw;


// construct & destruct
// --
CContactMesh::CContactMesh()
{
    ;
}
CContactMesh::~CContactMesh()
{
    // -> mgLevel別にdeleteしていく必要があるので,MeshのNodeと同様の処理でdeleteする.
    //
    vector<CContactNode*>::iterator itConNode;
    CContactNode *pConNode;
    for(itConNode=mvConNode.begin(); itConNode< mvConNode.end(); itConNode++){
        pConNode = *itConNode;
        if(mLevel==pConNode->getLevel()){
            delete pConNode;
            mvConNode.erase(itConNode);
            itConNode--;
        }
    };

    for_each(mvFace.begin(), mvFace.end(), DeleteObject());
    for_each(mvSlaveFace.begin(), mvSlaveFace.end(), DeleteObject());

    std::cout  << "~CContactMesh" << ", mgLevel= " << mLevel << std::endl;
}

// ContactNodeの追加, mapデータへID,Indexを代入
//
void CContactMesh::addConNode(CContactNode* pConNode, const uint& id)
{
    mvConNode.push_back(pConNode);
    
    mmConNodeID2Index[id]= mvConNode.size()-1;
}
void CContactMesh::addMasterConNode(CContactNode* pConNode, const uint& id)
{
    mvMasterConNode.push_back(pConNode);

    mmMasterConNodeID2Index[id]= mvMasterConNode.size()-1;
}
void CContactMesh::addSlaveConNode(CContactNode* pConNode, const uint& id)
{
    mvSlaveConNode.push_back(pConNode);

    mmSlaveConNodeID2Index[id]= mvSlaveConNode.size()-1;
}


// MasterFaceの追加, mapデータへID,Indexを代入
//
void CContactMesh::addMasterMeshID(const uint& id)
{
    mvMasterMeshID.push_back(id);

    uint index= mvMasterMeshID.size()-1;
    mmMasterMeshID2Index[id]= index;
}

// SlaveFaceの追加, mapデータへID,Indexを代入
//
void CContactMesh::addSlaveMeshID(const uint& id)
{
    mvSlaveMeshID.push_back(id);

    uint index= mvSlaveMeshID.size()-1;
    mmSlaveMeshID2Index[id]= index;
}


// ---
// Refine関連
// ---
// ConNodeに接続している(MasterFace,SlaveFace)のIDを集める.
//
void CContactMesh::setupAggSkinFace()
{
    CSkinFace  *pMFace;
    CSkinFace  *pSFace;
    CContactNode *pConNode;

    uint numOfMFace(mvFace.size()), numOfSFace(mvSlaveFace.size()), numOfConNode;
    uint iface,icnode;
    for(iface=0; iface< numOfMFace; iface++){
        pMFace= mvFace[iface];
        numOfConNode= pMFace->getNumOfNode();

        for(icnode=0; icnode< numOfConNode; icnode++){
            pConNode= pMFace->getNode(icnode);
            pConNode->setAggElemID(pMFace->getID());
            pConNode->setNeibElemVert(pMFace->getID(), icnode);
        };
    };
    for(iface=0; iface< numOfSFace; iface++){
        pSFace= mvSlaveFace[iface];
        numOfConNode= pSFace->getNumOfNode();

        for(icnode=0; icnode< numOfConNode; icnode++){
            pConNode= pSFace->getNode(icnode);
            pConNode->setAggElemID(pSFace->getID());
            pConNode->setNeibElemVert(pSFace->getID(), icnode);
        };
    };
}


// 辺に接続するSkinFaceのセットと
//    辺ノードの生成,設置を同時に行う.
// 
//  countID:新しい辺ノードのためのノードIDカウンター
//  pProgConMesh:プロロンゲーションする上位のContactMesh
// 
void CContactMesh::setupEdgeConNode(CContactMesh *pProgConMesh, uint& countID)
{
    PairConNode pairConNode;
    CSkinFace *pFace, *pNeibFace;
    CContactNode* pEdgeConNode;
    vdouble vCoord; vCoord.resize(3);
    uint numOfFace, numOfEdge;
    uint level;//EdgeConNodeのmgLevel
    uint numOfScalar,numOfDisp;//新ContactNodeのDOF
    uint iface,iedge;
    uint maslave;
    // 0:マスター面, 1:スレーブ面
    for(maslave=0; maslave< 2; maslave++){
        
        if(maslave==0) numOfFace= mvFace.size();     //マスター面ケース
        if(maslave==1) numOfFace= mvSlaveFace.size();//スレーブ面ケース
        
        for(iface=0; iface< numOfFace; iface++){
            if(maslave==0) pFace= mvFace[iface];     //マスター面ケース
            if(maslave==1) pFace= mvSlaveFace[iface];//スレーブ面ケース
            numOfEdge= pFace->getNumOfEdge();

            for(iedge=0; iedge< numOfEdge; iedge++){
                //未設置の辺であれば;
                //  ConNodeを生成しセット=> 隣接FaceにもConNodeをセット.
                //  隣接Faceをセット => 隣接Faceには自分Faceをセット.
                if(!pFace->isEdgeNodeMarking(iedge)){
                    pairConNode= pFace->getEdgePairNode(iedge);

                    uint numOfAggFaceA= pairConNode.first->getNumOfAggElem();
                    uint numOfAggFaceB= pairConNode.second->getNumOfAggElem();
                    uint faceIDa,faceIDb;
                    uint iagg,jagg;
                    bool bfind(false);//隣接Faceが見つかった場合は,Trueに設定
                    // 隣接するFaceを検索してセット,かつ,辺ConNodeを生成してセット
                    // ----
                    for(iagg=0; iagg< numOfAggFaceA; iagg++){
                        faceIDa= pairConNode.first->getAggElemID(iagg);

                        for(jagg=0; jagg< numOfAggFaceB; jagg++){
                            faceIDb= pairConNode.second->getAggElemID(jagg);
                            // 辺に接続するSkinFace(pFace自身は除外)
                            // --
                            if(faceIDa==faceIDb && faceIDa!=pFace->getID()){

//                                //debug
//                                cout << "CContactMesh::setupEdgeConNode, pFace ID= " << pFace->getID()
//                                     << ", 隣接Face ID= " << faceIDa
//                                     << ", Edge Index=  " << iedge
//                                     << endl;

                                //// Faceのインデックス: ID==Index として処理
                                if(maslave==0) pNeibFace= mvFace[faceIDa];     //マスター面ケース
                                if(maslave==1) pNeibFace= mvSlaveFace[faceIDa];//スレーブ面ケース
                                //// 辺ノード生成
                                pEdgeConNode= new CContactNode;// <<<<<<<<<- ContactNode生成
                                //------------------------
                                // 辺ノードへの座標,IDのセット
                                //------------------------
                                vCoord[0]= pairConNode.first->getX(); vCoord[0]+= pairConNode.second->getX(); vCoord[0] *= 0.5;
                                vCoord[1]= pairConNode.first->getY(); vCoord[1]+= pairConNode.second->getY(); vCoord[1] *= 0.5;
                                vCoord[2]= pairConNode.first->getZ(); vCoord[2]+= pairConNode.second->getZ(); vCoord[2] *= 0.5;
                                pEdgeConNode->setCoord(vCoord);
                                
                                level= pairConNode.first->getLevel() + 1;//頂点のConNodeのmgLevelより上段へ
                                pEdgeConNode->setLevel(level);
                                
                                pEdgeConNode->setID(countID);
                                countID++;////////// countIDの定義は,Factory,次の為にカウントアップ

                                //新ConNodeに,numOfScalar,numOfDispの設定
                                numOfScalar= pairConNode.first->getNumOfScalar();
                                numOfDisp= pairConNode.first->getNumOfDisp();
                                pEdgeConNode->resizeDisp(numOfDisp); pEdgeConNode->initDisp();
                                pEdgeConNode->resizeScalar(numOfScalar); pEdgeConNode->initScalar();
                                
                                //// 上位 接合メッシュにノードを追加
                                pProgConMesh->addConNode(pEdgeConNode, countID);//ContactMesh全体のConNode配列
                                if(maslave==0) pProgConMesh->addMasterConNode(pEdgeConNode, countID);//マスターConNode配列
                                if(maslave==1) pProgConMesh->addSlaveConNode(pEdgeConNode, countID);//スレーブConNode配列

                                // pFaceへの処理
                                pFace->setEdgeFace(pNeibFace, iedge);
                                pFace->setEdgeConNode(pEdgeConNode, iedge);
                                pFace->markingEdgeNode(iedge);
                                // pNeibFaceへの処理
                                pNeibFace->setEdgeFace(pFace, pairConNode);
                                pNeibFace->setEdgeConNode(pEdgeConNode, pairConNode);
                                pNeibFace->markingEdgeNode(pairConNode);

                                bfind=true;//1辺に隣接するFaceは一つなので(pFace自身は除外してる)
                                break;
                            }
                        };// jaggループ

                        if(bfind) break;//iaggのループをブレイク(既に隣接Faceを発見したので)

                    };// iaggループ

                    // 隣接していない辺にConNodeを生成
                    // ----
                    if(!bfind){
                        //// 辺ノード生成
                        pEdgeConNode= new CContactNode;// <<<<<<<<<- ContactNode生成
                        //------------------------
                        // 辺ノードへの座標,IDのセット
                        //------------------------
                        vCoord[0]= pairConNode.first->getX(); vCoord[0]+= pairConNode.second->getX(); vCoord[0] *= 0.5;
                        vCoord[1]= pairConNode.first->getY(); vCoord[1]+= pairConNode.second->getY(); vCoord[1] *= 0.5;
                        vCoord[2]= pairConNode.first->getZ(); vCoord[2]+= pairConNode.second->getZ(); vCoord[2] *= 0.5;
                        pEdgeConNode->setCoord(vCoord);

                        level= pairConNode.first->getLevel() + 1;//頂点のConNodeのmgLevelより上段へ
                        pEdgeConNode->setLevel(level);

                        pEdgeConNode->setID(countID);
                        countID++;////////// countIDの定義は,Factory,次の為にカウントアップ

                        //新ConNodeに,numOfScalar,numOfDispの設定
                        numOfScalar= pairConNode.first->getNumOfScalar();
                        numOfDisp= pairConNode.first->getNumOfDisp();
                        pEdgeConNode->resizeDisp(numOfDisp); pEdgeConNode->initDisp();
                        pEdgeConNode->resizeScalar(numOfScalar); pEdgeConNode->initScalar();

                        //// 上位 接合メッシュにノードを追加
                        pProgConMesh->addConNode(pEdgeConNode, countID);//ContactMesh全体のConNode配列
                        if(maslave==0) pProgConMesh->addMasterConNode(pEdgeConNode, countID);//マスターConNode配列
                        if(maslave==1) pProgConMesh->addSlaveConNode(pEdgeConNode, countID);//スレーブConNode配列

                        // pFaceへの処理
                        pFace->setEdgeConNode(pEdgeConNode, iedge);
                        pFace->markingEdgeNode(iedge);

//                        //debug
//                        cout << "CContactMesh::setupEdgeConNode, pFace ID= " << pFace->getID()
//                             << ", 隣接Faceなし, Edge Index=  " << iedge
//                             << endl;
                    }

                }// !isEdgeNodeMarking
            };// iedgeループ
        };// ifaceループ
    };//Master=0,Slave=1 ループ

    
    // Factory::Refineで処理..
    //
    // ConNodeへNodeIDのセット,FaceへElementIDのセット
    // --
    //  SkinFaceのmbSelfDomがTrueならば,Meshから辺ノードのID,要素IDを取得
    // 要素IDは,親要素からたどる.
}

// 面中心ノードの生成,設置
//  countID:新しく生成する面ノードのID
//  pProgConMesh:プロロンゲーションする上位のContactMesh
void CContactMesh::setupFaceConNode(CContactMesh *pProgConMesh, uint& countID)
{
    CSkinFace *pFace;
    CContactNode *pConNode, *pConFaceNode;
    uint numOfFace,numOfNode;
    uint numOfScalar,numOfDisp;//新ContactNodeのDOF
    uint iface,inode;
    vdouble vCoord; vCoord.resize(3);//面中心の座標

    uint maslave;
    for(maslave=0; maslave< 2; maslave++){
        if(maslave==0) numOfFace= mvFace.size();     //マスター面
        if(maslave==1) numOfFace= mvSlaveFace.size();//スレーブ面

        for(iface=0; iface< numOfFace; iface++){
            if(maslave==0) pFace= mvFace[iface];     //マスター面
            if(maslave==1) pFace= mvSlaveFace[iface];//スレーブ面

            numOfNode= pFace->getNumOfNode();
            vCoord[0]=0.0; vCoord[1]=0.0; vCoord[2]=0.0;
            for(inode=0; inode< numOfNode; inode++){
                pConNode= pFace->getNode(inode);
                vCoord[0]+= pConNode->getX(); vCoord[1]+= pConNode->getY(); vCoord[2]+= pConNode->getZ();
            };
            vCoord[0] /= (double)numOfNode;  vCoord[1] /= (double)numOfNode;  vCoord[2] /= (double)numOfNode;

            //// 面ノード生成
            pConFaceNode = new CContactNode;// <<<<<<<<<<<- 面中心のノード
            pConFaceNode->setCoord(vCoord);
            pConFaceNode->setID(countID);

            countID++;// 次の為にID番号をカウントアップ,countIDの定義は,Factory

            //新ConNodeに,numOfScalar,numOfDispの設定
            numOfScalar= pConNode->getNumOfScalar();//面構成ConNodeの最後のIndexのConNodeから自由度を取得
            numOfDisp= pConNode->getNumOfDisp();    //--
            pConFaceNode->resizeDisp(numOfDisp); pConFaceNode->initDisp();
            pConFaceNode->resizeScalar(numOfScalar); pConFaceNode->initScalar();

            //// 上位 接合メッシュにノードを追加
            pProgConMesh->addConNode(pConFaceNode,countID);//ContactMesh全体のConNode配列
            if(maslave==0) pProgConMesh->addMasterConNode(pConFaceNode, countID);//マスターConNode配列
            if(maslave==1) pProgConMesh->addSlaveConNode(pConFaceNode, countID);//スレーブConNode配列

            //// pFaceにpConFaceNodeをセット
            pFace->setFaceConNode(pConFaceNode);

        };//faceループ
    };//master,slave ループ
}


// 再分割した複数のFaceをセット
// --
// Master
void CContactMesh::addMasterFace(vector<CSkinFace*>& vface)
{
    CSkinFace* pFace;
    uint i, numOfFace(vface.size());
    for(i=0; i< numOfFace; i++){
        pFace= vface[i];
        mvFace.push_back(pFace);
    };
}
// Slave
void CContactMesh::addSlaveFace(vector<CSkinFace*>& vface)
{
    CSkinFace* pFace;
    uint i, numOfFace(vface.size());
    for(i=0; i< numOfFace; i++){
        pFace= vface[i];
        mvSlaveFace.push_back(pFace);
    };
}

// BBoxのみ利用バージョン(Octreeなし)
//  マスター面にスレーブ点を設置 -> Equationは,MasterFace関数
//
void CContactMesh::setupSPointOnMFace()
{
    uint numOfSlaveNode= mvSlaveConNode.size();
    CContactNode* pConNode;
    uint isnode;

    uint numOfMaster= mvFace.size();
    CSkinFace *pMFace;
    uint imface;
    for(imface=0; imface< numOfMaster; imface++){
        pMFace= mvFace[imface];

        pMFace->CalcNzNormalVector();//マスター面の面ベクトル(正規化)
        mBoundingBox.sizingOBB(pMFace);//マスター面に合わせたBBox

        //Octreeなしで,強引にBBoxで絞り込み
        // Octreeの場合は,全Slave点のループでは無くなる.
        for(isnode=0; isnode< numOfSlaveNode; isnode++){
            pConNode= mvSlaveConNode[isnode];

            if(mBoundingBox.judgeOBB(pConNode)) pMFace->addSlaveNode(pConNode);//スレーブ点の内外判定->内ならば,スレーブとして追加
        };//スレーブ点ループ
    };//マスター面ループ
}




