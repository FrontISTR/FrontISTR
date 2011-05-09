
#include "SkinFace.h"

//
// ContactMesh.cpp
//
//			2009.01.08
//			2009.01.08
//			k.Takeda
#include <vector>
#include "Vertex.h"
#include "ContactNode.h"
#include "OctreeKnot.h"
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
    std::cout  << "~CContactMesh, start" << std::endl;

//    // mgLevel別にdeleteしていく必要があるので,MeshのNodeと同様の処理でdeleteする.
//    //
//    vector<CContactNode*>::iterator itConNode;
//    CContactNode *pConNode;
//    for(itConNode=mvConNode.begin(); itConNode!=mvConNode.end(); itConNode++){
//        pConNode = *itConNode;
//        if(pConNode){
//            if(mLevel==pConNode->getLevel()){
//                std::cout << "aaaaa" << std::endl;
//                mvConNode.erase(itConNode);
//                std::cout << "bbbbb" << std::endl;
//                //itConNode--;
//                delete pConNode;
//                std::cout << "ccccc" << std::endl;
//            }
//        }
//    };
    
    // mgLevel別にdelete 2010.05.31 VC++仕様に合わせる
    //
    CContactNode *pConNode;
    uint i;
    for(i=0; i < mvConNode.size() ; i++){
        pConNode = mvConNode[i];
        if(mLevel==pConNode->getLevel()) delete pConNode;
            mvConNode.erase(mvConNode.begin()+i);
    };

    for_each(mvFace.begin(), mvFace.end(), DeleteObject());
    for_each(mvSlaveFace.begin(), mvSlaveFace.end(), DeleteObject());


    // mvKnot (Octreeオブジェクトの破棄)
    // 
    vector<COctreeKnot*> vKnot;
    uint numOfLayer= mvKnot.size();
    uint ilayer;
    // Layer==0 は moOctreeKnot なので除外
    for(ilayer=1; ilayer < numOfLayer; ilayer++){
        vKnot= mvKnot[ilayer];

        for_each(vKnot.begin(), vKnot.end(), DeleteObject());
    };


    std::cout  << "~CContactMesh, end " << ", mgLevel= " << mLevel << std::endl;
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

// 自身のConNodeを,progConMeshにセットアップ
//
void CContactMesh::setupCoarseConNode(CContactMesh* pProgConMesh)
{
    uint numOfNode;
    numOfNode= mvConNode.size();
    CContactNode *pConNode;

    uint inode, id;
    //全体
    for(inode=0; inode< numOfNode; inode++){
        pConNode= mvConNode[inode];
        pConNode->pushLevelMarking();//2010.05.27

        id= pConNode->getID();

        pProgConMesh->addConNode(pConNode, id);
    };

    //MasterConNode
    numOfNode= mvMasterConNode.size();
    for(inode=0; inode< numOfNode; inode++){
        pConNode= mvMasterConNode[inode];
        id= pConNode->getID();

        pProgConMesh->addMasterConNode(pConNode, id);
    };

    //SlaveConNode
    numOfNode= mvSlaveConNode.size();
    for(inode=0; inode< numOfNode; inode++){
        pConNode= mvSlaveConNode[inode];
        id= pConNode->getID();

        pProgConMesh->addSlaveConNode(pConNode, id);
    };
}

// ConNodeに接続している(MasterFace,SlaveFace)のIDを集める.
//
void CContactMesh::setupAggSkinFace()
{
    CSkinFace  *pMFace;
    CSkinFace  *pSFace;
    CContactNode *pConNode;
    
    uint numOfMFace(mvFace.size()), numOfSFace(mvSlaveFace.size()), numOfConNode;
    uint iface,icnode;

    // 第2段 以降のprolongationのときに,
    //  VertexのElemIndexに同じIDが入らないようにクリア. 10.03.16
    numOfConNode= mvMasterConNode.size();
    for(icnode=0; icnode< numOfConNode; icnode++){
        pConNode= mvMasterConNode[icnode];

        pConNode->clearAggElemID();
        pConNode->clearNeibElemVert();
    };
    numOfConNode= mvSlaveConNode.size();
    for(icnode=0; icnode< numOfConNode; icnode++){
        pConNode= mvSlaveConNode[icnode];

        pConNode->clearAggElemID();
        pConNode->clearNeibElemVert();
    };




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
void CContactMesh::setupEdgeConNode(CContactMesh *pProgConMesh)
{
    PairConNode pairConNode;
    CSkinFace *pFace, *pNeibFace;
    CContactNode* pEdgeConNode;
    vdouble vCoord; vCoord.resize(3);
    uint numOfFace, numOfEdge;
    uint level; //EdgeConNodeのmgLevel
    uint numOfScalar,numOfDisp;//新ContactNodeのDOF
    uint face_rank,rank;//新ConNodeのランク設定用
    uint iface,iedge;
    uint countID=mvConNode.size();//マスター,スレーブ関係なくID割り振り
    uint maslave;
    // 0:マスター面, 1:スレーブ面
    for(maslave=0; maslave< 2; maslave++){
        
        if(maslave==0) numOfFace= mvFace.size();     //マスター面ケース
        if(maslave==1) numOfFace= mvSlaveFace.size();//スレーブ面ケース
        
        for(iface=0; iface< numOfFace; iface++){
            if(maslave==0) pFace= mvFace[iface];     //マスター面ケース
            if(maslave==1) pFace= mvSlaveFace[iface];//スレーブ面ケース

            numOfEdge= pFace->getNumOfEdge();
            face_rank= pFace->getRank();

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

                                //// Faceのインデックス: ID==Index として処理
                                //if(maslave==0) pNeibFace= mvFace[faceIDa];     //マスター面ケース
                                //if(maslave==1) pNeibFace= mvSlaveFace[faceIDa];//スレーブ面ケース
                                // ----
                                // ID => Index を利用して隣接Faceを取得
                                // ----
                                uint index;
                                if(maslave==0){
                                    index= mmMasterFaceID2Index[faceIDa];
                                    pNeibFace= mvFace[index];
                                }
                                if(maslave==1){
                                    index= mmSlaveFaceID2Index[faceIDa];
                                    pNeibFace= mvSlaveFace[index];
                                }



                                //// 辺ノード生成
                                pEdgeConNode= new CContactNode;// <<<<<<<<<- ContactNode生成
                                //------------------------
                                // 辺ノードへの座標,IDのセット
                                //------------------------
                                vCoord[0]= pairConNode.first->getX(); vCoord[0]+= pairConNode.second->getX(); vCoord[0] *= 0.5;
                                vCoord[1]= pairConNode.first->getY(); vCoord[1]+= pairConNode.second->getY(); vCoord[1] *= 0.5;
                                vCoord[2]= pairConNode.first->getZ(); vCoord[2]+= pairConNode.second->getZ(); vCoord[2] *= 0.5;
                                pEdgeConNode->setCoord(vCoord);

                                
                                //level= pairConNode.first->getLevel() + 1;//頂点のConNodeのmgLevelより上段へ
                                level= mLevel+1;//2010.05.27
                                pEdgeConNode->setLevel(level);
                                pEdgeConNode->pushLevelMarking();//2010.05.27
                                pEdgeConNode->setID(countID);
                                countID++;////////// 次の為にカウントアップ
                                
                                //ConNodeが自身のMeshに含まれているConNodeならば,MeshIDをセット
                                if(pairConNode.first->getRank()==myRank){
                                    uint meshID= pairConNode.first->getMeshID();
                                    pEdgeConNode->setMeshID(meshID);
                                }

                                // ConNodeランク設定: 節点界面で考えてある.
                                if((pairConNode.first->getRank()!=face_rank) && (pairConNode.second->getRank()!=face_rank)){
                                    rank= pairConNode.first->getRank();
                                }else{
                                    rank= face_rank;
                                }
                                pEdgeConNode->setRank(rank);


                                //新ConNodeに,numOfScalar,numOfDispの設定
                                numOfScalar= pairConNode.first->getNumOfScalar();
                                numOfDisp= pairConNode.first->getNumOfDisp();
                                pEdgeConNode->resizeDisp(numOfDisp); pEdgeConNode->initDisp();
                                pEdgeConNode->resizeScalar(numOfScalar); pEdgeConNode->initScalar();
                                
                                //// 上位 接合メッシュにノードを追加
                                pProgConMesh->addConNode(pEdgeConNode, countID);                     //全体のConNode配列
                                if(maslave==0) pProgConMesh->addMasterConNode(pEdgeConNode, countID);//マスターConNode配列
                                if(maslave==1) pProgConMesh->addSlaveConNode(pEdgeConNode, countID); //スレーブConNode配列

                                // pFaceへの処理
                                pFace->setEdgeFace(pNeibFace, iedge);
                                pFace->setEdgeConNode(pEdgeConNode, iedge);
                                pFace->markingEdgeNode(iedge);
                                //pFace->setEdgeFace(pNeibFace, pairConNode);
                                //pFace->setEdgeConNode(pEdgeConNode, pairConNode);
                                //pFace->markingEdgeNode(pairConNode);

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

                        
                        //level= pairConNode.first->getLevel() + 1;//頂点のConNodeのmgLevelより上段へ
                        level= mLevel+1;//2010.05.27
                        pEdgeConNode->setLevel(level);
                        pEdgeConNode->pushLevelMarking();//2010.05.27
                        pEdgeConNode->setID(countID);
                        countID++;////////// 次の為にカウントアップ
                        
                        // ConNodeが自身のMeshに含まれているConNodeならば,MeshIDをセット
                        if(pairConNode.first->getRank()==myRank){
                            uint meshID= pairConNode.first->getMeshID();
                            pEdgeConNode->setMeshID(meshID);
                        }

                        // ConNodeランク設定: 節点界面で考えてある.
                        if((pairConNode.first->getRank()!=face_rank) && (pairConNode.second->getRank()!=face_rank)){
                            rank= pairConNode.first->getRank();
                        }else{
                            rank= face_rank;
                        }
                        pEdgeConNode->setRank(rank);

                        
                        //新ConNodeに,numOfScalar,numOfDispの設定
                        numOfScalar= pairConNode.first->getNumOfScalar();
                        numOfDisp= pairConNode.first->getNumOfDisp();
                        pEdgeConNode->resizeDisp(numOfDisp); pEdgeConNode->initDisp();
                        pEdgeConNode->resizeScalar(numOfScalar); pEdgeConNode->initScalar();

                        //// 上位 接合メッシュにノードを追加
                        pProgConMesh->addConNode(pEdgeConNode, countID);                     //全体のConNode配列
                        if(maslave==0) pProgConMesh->addMasterConNode(pEdgeConNode, countID);//マスターConNode配列
                        if(maslave==1) pProgConMesh->addSlaveConNode(pEdgeConNode, countID); //スレーブConNode配列

                        // pFaceへの処理
                        pFace->setEdgeConNode(pEdgeConNode, iedge);
                        pFace->markingEdgeNode(iedge);
                        //pFace->setEdgeConNode(pEdgeConNode, pairConNode);
                        //pFace->markingEdgeNode(pairConNode);
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
void CContactMesh::setupFaceConNode(CContactMesh *pProgConMesh)
{
    CSkinFace *pFace;
    CContactNode *pConNode, *pConFaceNode;
    uint numOfFace,numOfNode;
    uint numOfScalar,numOfDisp;//新ContactNodeのDOF
    uint iface,inode;
    uint countID= pProgConMesh->getNumOfConNode();//EdgeNodeも加算されたサイズを取得(マスター,スレーブ関係なくID割り振り)
    uint face_rank;
    vdouble vCoord; vCoord.resize(3);//面中心の座標

    uint maslave;
    for(maslave=0; maslave< 2; maslave++){
        if(maslave==0) numOfFace= mvFace.size();     //マスター面
        if(maslave==1) numOfFace= mvSlaveFace.size();//スレーブ面

        for(iface=0; iface< numOfFace; iface++){
            if(maslave==0) pFace= mvFace[iface];     //マスター面
            if(maslave==1) pFace= mvSlaveFace[iface];//スレーブ面

            face_rank= pFace->getRank();

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
            countID++;////////// 次の為にID番号をカウントアップ

            pConFaceNode->setRank(face_rank);//面中心の新ConNodeは,Faceのrankと同一(節点界面を採用)'10.02.28
            pConFaceNode->setLevel(mLevel+1);//一段上のLevel
            pConFaceNode->pushLevelMarking();//2010.05.27
            
            //ConNodeが自身のMeshに含まれているConNodeならば,MeshIDをセット
            if(pConNode->getRank()==myRank){
                uint meshID= pConNode->getMeshID();
                pConFaceNode->setMeshID(meshID);
            }


            //新ConNodeに,numOfScalar,numOfDispの設定
            numOfScalar= pConNode->getNumOfScalar();//面構成ConNodeの最後のIndexのConNodeから自由度を取得
            numOfDisp= pConNode->getNumOfDisp();    //--
            pConFaceNode->resizeDisp(numOfDisp); pConFaceNode->initDisp();
            pConFaceNode->resizeScalar(numOfScalar); pConFaceNode->initScalar();

            //// 上位 接合メッシュにノードを追加
            pProgConMesh->addConNode(pConFaceNode,countID);                      //全体のConNode配列
            if(maslave==0) pProgConMesh->addMasterConNode(pConFaceNode, countID);//マスターConNode配列
            if(maslave==1) pProgConMesh->addSlaveConNode(pConFaceNode, countID); //スレーブConNode配列

            //// pFaceにpConFaceNodeをセット
            pFace->setFaceConNode(pConFaceNode);

        };//faceループ
    };//master,slave ループ
}

//
//
void CContactMesh::addMasterFace(CSkinFace* pFace)
{
    mvFace.push_back(pFace);
    mmMasterFaceID2Index[pFace->getID()]= mvFace.size()-1;
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
        mmMasterFaceID2Index[pFace->getID()]= mvFace.size()-1;
    };
}
// Slave
void CContactMesh::addSlaveFace(CSkinFace* pSlaveFace)
{
    mvSlaveFace.push_back(pSlaveFace);
    mmSlaveFaceID2Index[pSlaveFace->getID()]= mvSlaveFace.size()-1;
}
void CContactMesh::addSlaveFace(vector<CSkinFace*>& vface)
{
    CSkinFace* pFace;
    uint i, numOfFace(vface.size());
    for(i=0; i< numOfFace; i++){
        pFace= vface[i];
        mvSlaveFace.push_back(pFace);
        mmSlaveFaceID2Index[pFace->getID()]= mvSlaveFace.size()-1;
    };
}


//  マスター面にスレーブ点を設置 -> Equationは,MasterFace関数
//
void CContactMesh::setupSPointOnMFace()
{

    // Octreeによる検索
    //
    uint maxLayer= mvKnot.size()-1;
    uint numOfMaster= mvFace.size();
    CSkinFace *pMFace;
    uint imface;
    
    ////debug
    //ofstream ofsd;
    //ofsd.open("setupSPointOnMFace.out", ios::out | ios::app);//ファイルの後ろに追加
    //ofsd.open("setupSPointOnMFace.out", ios::out);

    for(imface=0; imface< numOfMaster; imface++){

        pMFace= mvFace[imface];
        
        //マスター面すべてが,同じOctreeKnotに入るまで,
        //  Octreeのレイヤーを上げて,Octreeレイヤーを決定
        //
        CContactNode *pVertNode;
        uint numOfVert= pMFace->getNumOfNode();
        uint sum_knotID, knotID;
        uint finalLayer(0);
        uint iLayer, ivert;
        for(iLayer=maxLayer; iLayer >= 0; iLayer--){

            sum_knotID =0;
            for(ivert=0; ivert < numOfVert; ivert++){
                pVertNode= pMFace->getNode(ivert);
                
                if(ivert==0) knotID= pVertNode->getOctreeID(iLayer);
                sum_knotID += pVertNode->getOctreeID(iLayer);
            };

            if(sum_knotID == knotID*numOfVert){//頂点の所属KnotIDが一致したレイヤー番号を取得
                finalLayer= iLayer;
                break;
            }
        };
        //マスター面が全て入るOctreeKnot
        COctreeKnot *pKnot= mvKnot[finalLayer][knotID];
        uint numOfSlaveNode= pKnot->getNumOfSlaveNode();
        
        
        pMFace->CalcNzNormalVector();  //マスター面の面ベクトル(正規化)
        mBoundingBox.sizingOBB(pMFace);//BBox(マスター面サイズ)

        CContactNode* pConNode;
        uint isnode;
        for(isnode=0; isnode< numOfSlaveNode; isnode++){

            pConNode= pKnot->getSlaveNode(isnode);

            //BBoxに入るスレーブ点か判定 => マスター面の中にあるか判定の上でスレーブとして追加
            if(mBoundingBox.judgeOBB(pConNode)){
                //debug
                cout << "ContactMesh::setupSPointOnMFace, pConNode id= " << pConNode->getID()
                     << ", Node id= " << pConNode->getNodeID() << endl;

                pMFace->addSlaveNode(pConNode);
            }
        };
    };
    ////debug
    //ofsd.close();


/*
    // 総当たりによる検索
    //
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

        for(isnode=0; isnode< numOfSlaveNode; isnode++){
            pConNode= mvSlaveConNode[isnode];

            //スレーブ点の内外判定->内ならば,スレーブとして追加
            if(mBoundingBox.judgeOBB(pConNode)){
                pMFace->addSlaveNode(pConNode);

                //debug
                if(pConNode->getID()==22)
                    cout << "conID==22: pMFace ID= " << pMFace->getID() << endl;
            }
        };
    };
*/

}


// スレーブ点に対する,マスター面構成ノードのCoef 計算
//
void CContactMesh::setupMPC_Coef()
{
    uint numOfMaster= mvFace.size();
    CSkinFace *pMFace;
    uint imface;

    uint numOfSlave;
    uint islave;

    for(imface=0; imface < numOfMaster; imface++){

        pMFace= mvFace[imface];
        
        numOfSlave= pMFace->getNumOfSlaveNode();
        //debug
        cout << "ContactMesh::setupMPC_Coef, numOfSlaveNode= " << numOfSlave << endl;

        for(islave=0; islave < numOfSlave; islave++){
            pMFace->CalcSlave(islave, MPCValueType::Displacement);
        };
    };
}

//// スレーブ点を所有しているMasterFaceのIDを返す.
////   => MasterFaceを構成しているConNodeにCoefが入っている.
////
//int CContactMesh::getMFaceID_at_Slave(const uint& islave)
//{
//    CContactNode *pConNode;
//
//    pConNode= mvSlaveConNode[islave];
//
//    if(pConNode->have_MasterFaceID()){
//        return (int)pConNode->getMasterFaceID();
//    }else{
//        return -1;
//    }
//}

// 
// id からマスター面を返す.
//
CSkinFace* CContactMesh::getMasterFace_ID(const uint& id)
{
    uint index= mmMasterFaceID2Index[id];

    return mvFace[index];
}


// 八分木全体の生成 (指定Layer数の八分木を生成)
//
void CContactMesh::generateOctree(const uint& maxLayer)
{
    //親八分木のサイズ
    double minX,maxX;
    double minY,maxY;
    double minZ,maxZ;
    
    CContactNode *pConNode;
    uint numOfConNode= mvConNode.size();
    uint inode;
    for(inode=0; inode< numOfConNode; inode++){

        pConNode= mvConNode[inode];
        pConNode->resizeOctreeID(maxLayer+1);

        if(inode==0){
            //初期値
            minX= pConNode->getX(); maxX= pConNode->getX();
            minY= pConNode->getY(); maxY= pConNode->getY();
            minZ= pConNode->getZ(); maxZ= pConNode->getZ();
        }else{
            if(minX > pConNode->getX()) minX= pConNode->getX();
            if(minY > pConNode->getY()) minY= pConNode->getY();
            if(minZ > pConNode->getZ()) minZ= pConNode->getZ();

            if(maxX < pConNode->getX()) maxX= pConNode->getX();
            if(maxY < pConNode->getY()) maxY= pConNode->getY();
            if(maxZ < pConNode->getZ()) maxZ= pConNode->getZ();
        }
    };
    //親八分木のサイズ
    // 1.座標軸に平行なサイズの場合オフセット
    // 2.その他:サイズをα分だけ拡大
    double difX= minX - maxX;
    double difY= minY - maxY;
    double difZ= minZ - maxZ;
    double Alpha(0.1);
    if(abs(difX) <= 1.0e-4) maxX += Alpha;
    if(abs(difY) <= 1.0e-4) maxY += Alpha;
    if(abs(difZ) <= 1.0e-4) maxZ += Alpha;
    moOctreeKnot.setX(minX-Alpha, maxX+Alpha);
    moOctreeKnot.setY(minY-Alpha, maxY+Alpha);
    moOctreeKnot.setZ(minZ-Alpha, maxZ+Alpha);
    //親八分木のレイヤー,箱ID番号
    moOctreeKnot.setLayerID(0);
    moOctreeKnot.setID(0);



    uint numOfMasterNode, numOfSlaveNode;
    numOfMasterNode= mvMasterConNode.size();
    numOfSlaveNode= mvSlaveConNode.size();

    moOctreeKnot.reserveMasterNode(numOfMasterNode);
    moOctreeKnot.reserveSlaveNode(numOfSlaveNode);
    
    //親八分木のマスター面の頂点,スレーブ点のセット
    //
    for(inode=0; inode < numOfMasterNode; inode++){
        pConNode= mvMasterConNode[inode];
        moOctreeKnot.addMasterNode(pConNode);
    };
    for(inode=0; inode < numOfSlaveNode; inode++){
        pConNode= mvSlaveConNode[inode];
        moOctreeKnot.addSlaveNode(pConNode);
    };
    //八分木のマスター,スレーブ点にレイヤーとOctreeIDをセット
    moOctreeKnot.setItemProp();



    //レイヤー別のOctreeKnotの生成
    COctreeKnot *pPrevKnot,*pNextKnot;
    vector<COctreeKnot*> nextKnot;
    vector<COctreeKnot*> prevKnot;
    
    mvKnot.resize(maxLayer+1);         //"maxLayer+1" 分のLayer数
    mvKnot[0].push_back(&moOctreeKnot);//Layer==0のKnotは,親Knotの1個

    uint lastPos;
    uint ilayer, prev_pos, child_pos;
    // "maxLayer+1" ぶんのレイヤー数
    //
    for(ilayer=0; ilayer < maxLayer; ilayer++){
        
        uint knotID(0); //nextKnotのID

        prevKnot= mvKnot[ilayer];      //プレブKnot
        lastPos= mvKnot[ilayer].size();//prevKnot ラストカウント
        
        nextKnot.clear();
        for(prev_pos=0; prev_pos < lastPos; prev_pos++){
            
            pPrevKnot= prevKnot[prev_pos];
            pPrevKnot->createChildKnot();
            pPrevKnot->distItem();

            for(child_pos=0; child_pos < 8; child_pos++){

                pNextKnot= pPrevKnot->getChildKnot(child_pos);
                pNextKnot->setLayerID(ilayer+1);
                pNextKnot->setID(knotID);
                knotID++;

                pNextKnot->setItemProp();
                
                nextKnot.push_back(pNextKnot);
            };
        };
        mvKnot[ilayer+1]= nextKnot;
    };
    
}


// SkinFaceの辺ノードvectorのメモリー解放
//
void CContactMesh::deleteProgData()
{
    uint nNumOfFace;
    uint iFace;

    nNumOfFace = mvFace.size();
    for(iFace=0; iFace < nNumOfFace; iFace++) mvFace[iFace]->deleteProgData();

    nNumOfFace = mvSlaveFace.size();
    for(iFace=0; iFace < nNumOfFace; iFace++) mvSlaveFace[iFace]->deleteProgData();
}







