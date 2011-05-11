//
//  BoundaryVolumeMesh.cpp
//
//
//
//
//              2010.04.12
//              k.Takeda
#include <vector>

#include "BndVertex.h"
#include "BoundaryMesh.h"

#include "BoundaryVolume.h"

#include "BoundaryParts.h"
#include "BoundaryVolumeMesh.h"
using namespace pmw;


CBoundaryVolumeMesh::CBoundaryVolumeMesh()
{
    ;
}
CBoundaryVolumeMesh::~CBoundaryVolumeMesh()
{
    for_each(mvBVolume.begin(), mvBVolume.end(), DeleteObject());
}


void CBoundaryVolumeMesh::resizeVolume(const uiint& res_size)
{
    mvBVolume.resize(res_size);
}


void CBoundaryVolumeMesh::setBVolume(const uiint& index, CBoundaryVolume* pBVolume)
{
    uiint id;
    id= pBVolume->getID();

    mmBVolumeID2Index[id]= index;

    mvBVolume[index]= pBVolume;
}


void CBoundaryVolumeMesh::addBVolume(CBoundaryVolume* pBVolume)
{
    mvBVolume.push_back(pBVolume);

    uiint id;
    id= pBVolume->getID();

    mmBVolumeID2Index[id]= mvBVolume.size()-1;
}



// 頂点のVolID集合の領域確保
// ----
void CBoundaryVolumeMesh::resizeAggVol()
{
    uiint numOfAgg= mvBNode.size();
    
    mvAggregateVol.resize(numOfAgg);
}

// 頂点へVolumeIDを集合
// ----
void CBoundaryVolumeMesh::setupAggVol()
{
    uiint numOfBNode= mvBNode.size();
    CBoundaryNode  *pBNode;
    uiint ibnode;
    //下位グリッドのAggデータをクリア
    //----
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= mvBNode[ibnode];

        pBNode->clearAggElemID();
        pBNode->clearNeibElemVert();
    };
    
    uiint numOfVol= mvBVolume.size();
    CBoundaryVolume *pBVol;
    uiint ivol;
    //BNodeの周囲Volを集める
    //----
    for(ivol=0; ivol < numOfVol; ivol++){
        pBVol= mvBVolume[ivol];
        
        //uint numOfVert= pBVol->getNumOfBNode();
        uiint numOfVert= pBVol->getNumOfVert();
        uiint ivert;
        for(ivert=0; ivert < numOfVert; ivert++){
            pBNode= pBVol->getBNode(ivert);
            pBNode->setAggElemID(pBVol->getID());
        };
    };
    
    //Aggregateに移す.
    //----
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= mvBNode[ibnode];

        uiint numOfAgg= pBNode->getNumOfAggElem();
        uiint iAgg;

        for(iAgg=0; iAgg < numOfAgg; iAgg++){
            mvAggregateVol[ibnode].push_back(pBNode->getAggElemID(iAgg));
        };
    };
}


// EdgeBNode生成
//
void CBoundaryVolumeMesh::GeneEdgeBNode()
{
    uiint numOfVol= mvBVolume.size();
    uiint ivol;
    CBoundaryVolume *pBVol;

    CElement *pElem;      // 各BVolumeが所有するElement
    CNode *pNode0,*pNode1;// 辺の両端のNode(MeshのNode)
    CNode *pEdgeNode;     // 辺中央のNode

    uiint countID= mvBNode.size();//EdgeBNodeのID初期値

    for(ivol=0; ivol < numOfVol; ivol++){
        pBVol= mvBVolume[ivol];

        uiint numOfEdge= pBVol->getNumOfEdge();
        uiint iedge;
        PairBNode pairBNode;
        for(iedge=0; iedge < numOfEdge; iedge++){

            if(!pBVol->isMarkingEdge(iedge)){
                bool bfind(false);//辺に隣接Volが存在するか否かの真偽

                pairBNode= pBVol->getPairBNode(iedge);/// <<---
                
                // 辺中央のNodeを取得 => のちにEdgeBNodeにセット
                // ----
                pElem = pBVol->getElement();
                pNode0 = pairBNode.first->getNode();
                pNode1 = pairBNode.second->getNode();
                uiint nElemEdge= pElem->getEdgeIndex(pNode0, pNode1);


                uiint numOfiAgg= pairBNode.first->getNumOfAggElem();
                uiint numOfjAgg= pairBNode.second->getNumOfAggElem();
                uiint iAgg,jAgg, iVolID, jVolID;
                
                for(iAgg=0; iAgg < numOfiAgg; iAgg++){
                    iVolID= pairBNode.first->getAggElemID(iAgg);

                    for(jAgg=0; jAgg < numOfjAgg; jAgg++){
                        jVolID= pairBNode.second->getAggElemID(jAgg);

                        if(iVolID==jVolID){
                            if(iVolID != pBVol->getID()){

                                bfind= true;

                                // 生成
                                CBoundaryNode *pEdgeBNode= new CBoundaryNode;// <<<<<<<<<<<<<< new
                                
                                pEdgeBNode->setID(countID);
                                countID++;

                                // VolumeMeshにEdgeBNodeセット
                                mvBEdgeBNode.push_back(pEdgeBNode);

                                // pBVolの"辺"にEdgeBNodeをセット
                                pBVol->setEdgeNeibVol(iedge, iVolID);
                                pBVol->markingEdge(iedge);
                                pBVol->setEdgeBNode(iedge, pEdgeBNode);
                                

                                // 隣接するVolumeを取得
                                uiint neibIndex= mmBVolumeID2Index[iVolID];
                                CBoundaryVolume *pNeibBVol= mvBVolume[neibIndex];
                                
                                //cout << "BoundaryVolumeMesh::GeneEdgeBNode, ---- C 0" << endl;
                                uiint jedge= pNeibBVol->getEdgeID(pairBNode);
                                //cout << "BoundaryVolumeMesh::GeneEdgeBNode, ---- C 1" << endl;

                                // pNeibBVol(隣接Vol)の"辺"にEdgeBNodeをセット
                                pNeibBVol->setEdgeNeibVol(jedge, pBVol->getID());
                                pNeibBVol->markingEdge(jedge);
                                pNeibBVol->setEdgeBNode(jedge, pEdgeBNode);
                                

                                //辺中央のNodeを取得してEdgeBNodeにセット:for iAggの外に出した
                                //----
                                //pElem = pBVol->getElement();
                                //pNode0 = pairBNode.first->getNode();
                                //pNode1 = pairBNode.second->getNode();
                                //uint nEdgeIndex = pElem->getEdgeIndex(pNode0, pNode1);

                                pEdgeNode = pElem->getEdgeInterNode(nElemEdge);
                                pEdgeBNode->setNode(pEdgeNode);

                                //2次要素の場合 => 自身のmvBNodeに追加
                                //----
                                if(pBVol->getOrder()==ElementOrder::Second){
                                    uiint renum = mvBNode.size();
                                    mvBNode.resize(renum+1);
                                    mvBNode[renum]=pEdgeBNode;

                                    mmBNodeID2Index[pEdgeBNode->getID()] = renum;

                                    //辺BNodeに対応するAggIDをセット
                                    mvAggregateVol.resize(renum+1);
                                    mvAggregateVol[renum].push_back(pBVol->getID());
                                    mvAggregateVol[renum].push_back(pNeibBVol->getID());
                                    
                                    pBVol->replaceEdgeBNode(iedge);    //2次要素の場合: 要素内のEdgeBNode -> mvBNodeに移設
                                    pNeibBVol->replaceEdgeBNode(jedge);//2次要素の場合: 要素内のEdgeBNode -> mvBNodeに移設

                                    // EdgeBNodeの属性設定(2次要素)
                                    pEdgeBNode->setMGLevel(mMGLevel);
                                    pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel + 1);
                                }else{
                                    // EdgeBNodeの属性設定(1次要素)
                                    pEdgeBNode->setMGLevel(mMGLevel+1);
                                    pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel);

                                    mnEdgeNodeCount++;//progBMeshのmvBNodeサイズカウント用:this->refineで利用
                                }
                            }
                        }

                        if(bfind) break;
                    };// jAgg loop

                    if(bfind) break;
                };// iAgg loop

                // 隣接Volが存在しない場合
                if(!bfind){
                    // 生成
                    CBoundaryNode *pEdgeBNode= new CBoundaryNode;// <<<<<<<<<<<<<<< new
                    
                    pEdgeBNode->setID(countID);
                    countID++;

                    // VolumeMeshにEdgeBNodeセット
                    mvBEdgeBNode.push_back(pEdgeBNode);

                    // pBVolの"辺"にEdgeBNodeをセット
                    pBVol->markingEdge(iedge);
                    pBVol->setEdgeBNode(iedge, pEdgeBNode);


                    //辺中央のNodeを取得してEdgeBNodeにセット:for iAggの外に出した
                    //----
                    //pElem = pBVol->getElement();
                    //pNode0 = pairBNode.first->getNode();
                    //pNode1 = pairBNode.second->getNode();
                    //uint nEdgeIndex = pElem->getEdgeIndex(pNode0, pNode1);
                    
                    pEdgeNode = pElem->getEdgeInterNode(nElemEdge);
                    pEdgeBNode->setNode(pEdgeNode);


                    //2次要素の場合 => 自身のmvBNodeに追加
                    //----
                    if(pBVol->getOrder()==ElementOrder::Second){
                        uiint renum = mvBNode.size();
                        mvBNode.resize(renum+1);
                        mvBNode[renum]=pEdgeBNode;

                        mmBNodeID2Index[pEdgeBNode->getID()] = renum;

                        //辺BNodeに対応するAggIDをセット
                        mvAggregateVol.resize(renum+1);
                        mvAggregateVol[renum].push_back(pBVol->getID());

                        pBVol->replaceEdgeBNode(iedge);//2次要素の場合: 要素内のEdgeBNode -> mvBNodeに移設

                        // EdgeBNodeの属性設定(2次要素)
                        pEdgeBNode->setMGLevel(mMGLevel);
                        pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel + 1);
                    }else{
                        // EdgeBNodeの属性設定(1次要素)
                        pEdgeBNode->setMGLevel(mMGLevel+1);
                        pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel);

                        mnEdgeNodeCount++;//progBMeshのmvBNodeサイズカウント用:this->refineで利用
                    }
                }
            }// if(isMarkingEdge)
        };// iedge loop

    };// ivol loop
}

// FaceBNode生成
//
void CBoundaryVolumeMesh::GeneFaceBNode()
{
    uiint numOfVol= mvBVolume.size();
    CBoundaryVolume *pBVol;
    CBoundaryVolume *pNeibBVol;

    CElement *pElem;//各BVolumeが所有するElement
    vector<CNode*> vNode; vNode.resize(3);//面構成Node(MeshのNode)
    CNode *pFaceNode;//面中央のNode

    uiint countID= mvBNode.size() + mvBEdgeBNode.size();//FaceBNodeのID初期値
    uiint ivol;
    
    // 体積
    for(ivol=0; ivol < numOfVol; ivol++){
        
        pBVol= mvBVolume[ivol];
        
        uiint numOfFace= pBVol->getNumOfFace();
        uiint iface;
        vector<CBoundaryNode*> vBNode;

        // 面構成BNode
        for(iface=0; iface < numOfFace; iface++){

            if(!pBVol->isMarkingFace(iface)){

                vBNode= pBVol->getFaceCnvNodes(iface);

                uiint iVolID, jVolID, kVolID;

                uiint iAgg, jAgg, kAgg;
                uiint numOfiAgg= vBNode[0]->getNumOfAggElem();
                uiint numOfjAgg= vBNode[1]->getNumOfAggElem();
                uiint numOfkAgg= vBNode[2]->getNumOfAggElem();

                bool bfind(false);//隣接Volの存在

                
                // 隣接Volの検索
                for(iAgg=0; iAgg < numOfiAgg; iAgg++){
                    iVolID= vBNode[0]->getAggElemID(iAgg);

                    for(jAgg=0; jAgg < numOfjAgg; jAgg++){
                        jVolID= vBNode[1]->getAggElemID(jAgg);

                        for(kAgg=0; kAgg < numOfkAgg; kAgg++){
                            kVolID= vBNode[2]->getAggElemID(kAgg);

                            if(iVolID==jVolID && jVolID==kVolID && kVolID!=pBVol->getID()){

                                bfind= true;

                                CBoundaryNode *pFaceBNode= new CBoundaryNode;// <<<<<<<<<<<<<<<  new

                                // FaceBNodeの属性設定
                                pFaceBNode->setMGLevel(mMGLevel+1);
                                pFaceBNode->resizeValue(mMaxMGLevel-mMGLevel);
                                pFaceBNode->setID(countID);
                                countID++;
                                
                                mvBFaceBNode.push_back(pFaceBNode);

                                // Vol自身
                                pBVol->setFaceNeibVol(iface, iVolID);
                                pBVol->markingFace(iface);
                                pBVol->setFaceBNode(iface, pFaceBNode);


                                // 隣接Vol
                                uiint neibIndex= mmBVolumeID2Index[iVolID];
                                pNeibBVol= mvBVolume[neibIndex];
                                uiint neibFace= pNeibBVol->getFaceID(vBNode);

                                pNeibBVol->setFaceNeibVol(neibFace, pBVol->getID());
                                pNeibBVol->markingFace(neibFace);
                                pNeibBVol->setFaceBNode(neibFace, pFaceBNode);

                                
                                // Nodeのセット
                                // Elementから面Nodeを取得して、FaceBNodeにセット
                                // ----
                                pElem = pBVol->getElement();
                                uiint ivert;
                                for(ivert=0; ivert < 3; ivert++){
                                    vNode[ivert] = vBNode[ivert]->getNode();
                                };
                                uiint nFaceIndex = pElem->getFaceIndex(vNode[0],vNode[1],vNode[2]);
                                pFaceNode = pElem->getFaceNode(nFaceIndex);

                                pFaceBNode->setNode(pFaceNode);
                            }
                            if(bfind) break;
                        };// kAgg
                        if(bfind) break;
                    };// jAgg
                    if(bfind) break;
                };// iAgg


                // 隣接Volが無い面
                if(!bfind){
                    CBoundaryNode *pFaceBNode= new CBoundaryNode;// <<<<<<<<<<<<<<<<<< new

                    // FaceBNodeの属性設定
                    pFaceBNode->setMGLevel(mMGLevel+1);
                    pFaceBNode->resizeValue(mMaxMGLevel-mMGLevel);
                    pFaceBNode->setID(countID);
                    countID++;

                    mvBFaceBNode.push_back(pFaceBNode);

                    // Vol自身
                    pBVol->markingFace(iface);
                    pBVol->setFaceBNode(iface, pFaceBNode);

                    // Nodeのセット
                    // Elementから面Nodeを取得して、FaceBNodeにセット
                    // ----
                    pElem = pBVol->getElement();
                    uiint ivert;
                    for(ivert=0; ivert < 3; ivert++){
                        vNode[ivert] = vBNode[ivert]->getNode();
                    };
                    uiint nFaceIndex = pElem->getFaceIndex(vNode[0],vNode[1],vNode[2]);
                    pFaceNode = pElem->getFaceNode(nFaceIndex);

                    pFaceBNode->setNode(pFaceNode);
                }

            }// if(isMarkingFace)

        };// iface 
    };// ivol
}

// VolumeBNode生成
//
void CBoundaryVolumeMesh::GeneVolBNode()
{
    uiint numOfVol= mvBVolume.size();
    CBoundaryVolume *pBVol;
    CElement *pElem;//各BVoluemが所有するElement
    CNode *pNode;//Element中央のNode

    uiint countID= mvBNode.size() + mvBEdgeBNode.size() + mvBFaceBNode.size();//VolumeBNodeのID初期値
    uiint ivol;

    for(ivol=0; ivol < numOfVol; ivol++){
        pBVol= mvBVolume[ivol];

        CBoundaryNode *pBVolBNode= new CBoundaryNode;// <<<<<<<<<<<<<<<<<<< new

        pBVolBNode->setMGLevel(mMGLevel+1);
        pBVolBNode->resizeValue(mMaxMGLevel-mMGLevel);
        pBVolBNode->setID(countID);
        countID++;

        pBVol->setVolBNode(pBVolBNode);

        //Nodeのセット
        pElem = pBVol->getElement();
        pNode = pElem->getVolumeNode();
        pBVolBNode->setNode(pNode);
        

        mvBVolBNode.push_back(pBVolBNode);
    };
}


// ----
// 境界要素の再分割 => progBMeshにセット
// ----
void CBoundaryVolumeMesh::refine(CBoundaryVolumeMesh* pProgBVolMesh)
{
    CBoundaryVolume *pBVol;
    uiint ivol, numOfVol= mvBVolume.size();
    vector<CBoundaryVolume*> vProgVol;// pBVolから生成される子BVolume
    //----
    //BVolumeの分割 => progBVolMeshにセット
    //----
    uiint countID(0);//新ID
    for(ivol=0; ivol < numOfVol; ivol++){
        pBVol= mvBVolume[ivol];

        pBVol->refine(countID, mvDOF);//  Refine, 内部でcountID++

        vProgVol.clear();
        vProgVol= pBVol->getProgParts();

        uiint i;
        for(i=0; i < vProgVol.size(); i++){
            pProgBVolMesh->addBVolume(vProgVol[i]);
        };
    };

    //----
    //BNode(頂点,辺,面) => progBVolMeshにセット
    //----
    uiint numOfBNode    = mvBNode.size();
    //uint numOfEdgeBNode= mvBEdgeBNode.size();
    uiint numOfFaceBNode= mvBFaceBNode.size();
    uiint numOfVolBNode = mvBVolBNode.size();
    uiint numOfProgBNode= numOfBNode + mnEdgeNodeCount + numOfFaceBNode + numOfVolBNode;

    pProgBVolMesh->resizeBNode(numOfProgBNode);
    
    uiint ibnode;
    // 頂点BNode
    uiint init = 0;
    uiint end  = numOfBNode;
    
    for(ibnode=init; ibnode < end; ibnode++){
        pProgBVolMesh->setBNode(ibnode, mvBNode[ibnode]);
    };
    
    // 辺BNode
    init = numOfBNode;
    end  = numOfBNode + mnEdgeNodeCount;
    
    for(ibnode=init; ibnode < end; ibnode++){
        pProgBVolMesh->setBNode(ibnode, mvBEdgeBNode[ibnode-init]);
    };
    
    // 面BNode
    init = end;
    end += numOfFaceBNode;

    for(ibnode=init; ibnode < end; ibnode++){
        pProgBVolMesh->setBNode(ibnode, mvBFaceBNode[ibnode-init]);
    };

    // 体積BNode
    init = end;
    end += numOfVolBNode;

    for(ibnode=init; ibnode < end; ibnode++){
        pProgBVolMesh->setBNode(ibnode, mvBVolBNode[ibnode-init]);
    }
}


// ノイマン型境界値のBNodeへの再配分
//
void CBoundaryVolumeMesh::distNeumannValue()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    // ノイマン条件が設定されていなければ,Errorを返す.
    // --
    if(mnBndType != BoundaryType::Neumann){
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryVolumeMesh::distNeumannValue");
        return;
    }

    //形状関数
    CShapeHexa  *pShHexa  = CShapeHexa::Instance();
    CShapeTetra *pShTetra = CShapeTetra::Instance();
    CShapePrism *pShPrism = CShapePrism::Instance();

    uiint inode, numOfBNode=mvBNode.size();
    CBoundaryNode *pBNode;
    uiint dof, idof;
    //節点力-ゼロクリア
    for(inode=0; inode < numOfBNode; inode++){
        pBNode= mvBNode[inode];

        //uint numOfDOF= pBNode->getNumOfDOF();
        uiint numOfDOF = getNumOfDOF();

        for(idof=0; idof < numOfDOF; idof++){
            //dof= pBNode->getDOF(idof);
            dof = getDOF(idof);
            pBNode->initValue(dof, mMGLevel);//mgLevel別-自由度別の値をゼロクリア
        };
    };

    //形状関数の積分値による配分
    double integVal, nodalVal, entVal;
    uiint ivol, numOfVol=mvBVolume.size();
    uiint numOfDOF;
    for(ivol=0; ivol < numOfVol; ivol++){
        CBoundaryVolume *pBVol = mvBVolume[ivol];

        numOfDOF = getNumOfDOF();

        ////debug
        //cout << "CBoundaryVolumeMesh::distNeumannValue, mgLevel=" << mMGLevel << ", numOfDOF=" << numOfDOF << endl;

        // 自由度別の外力配分
        // ----
        for(idof=0; idof < numOfDOF; idof++){

            dof = getDOF(idof);//自由度番号
            entVal = pBVol->getBndValue(dof);//体積の所有する境界値

            ////debug
            //cout << "CBoundaryVolumeMesh::distNeumannValue, mgLevel=" << mMGLevel << ", dof=" << dof << ", entVal=" << entVal << endl;

            uiint ivert;
            switch(pBVol->getElemType()){
                case(ElementType::Hexa):
                    for(ivert=0; ivert < 8; ivert++){
                        integVal= pShHexa->getIntegralValue8(ivert);
                        nodalVal= entVal * integVal;// 等価節点力計算
                        
                        ////debug
                        //cout << "CBoundaryVolumeMesh::distNeumannValue, mgLevel=" << mMGLevel << ", dof=" << dof << ", nodalVal=" << nodalVal << endl;
                        
                        pBNode= pBVol->getBNode(ivert);// 境界Node
                        pBNode->addValue(dof, mMGLevel, nodalVal);// 境界値として節点"外力に加算"
                    };
                    break;
                case(ElementType::Hexa2):
                    for(ivert=0; ivert < 20; ivert++){
                        integVal= pShHexa->getIntegralValue20(ivert);
                        nodalVal= entVal * integVal;// 等価節点力計算

                        pBNode= pBVol->getBNode(ivert);// 境界Node
                        pBNode->addValue(dof, mMGLevel, nodalVal);// 境界値として節点"外力に加算"
                    };
                    break;
                case(ElementType::Tetra):
                    for(ivert=0; ivert < 4; ivert++){
                        integVal= pShTetra->getIntegValue4(ivert);
                        nodalVal= entVal * integVal;// 等価節点力
                        
                        pBNode= pBVol->getBNode(ivert);
                        pBNode->addValue(dof, mMGLevel, nodalVal);//外力に加算
                    }
                    break;
                case(ElementType::Tetra2):
                    for(ivert=0; ivert < 10; ivert++){
                        integVal= pShTetra->getIntegValue10(ivert);
                        nodalVal= entVal * integVal;//等価節点力

                        pBNode= pBVol->getBNode(ivert);
                        pBNode->addValue(dof, mMGLevel, nodalVal);//外力に加算
                    }
                    break;
                case(ElementType::Prism):
                    for(ivert=0; ivert < 6; ivert++){
                        integVal= pShPrism->getIntegValue6(ivert);
                        nodalVal= entVal * integVal;//等価節点力

                        pBNode= pBVol->getBNode(ivert);
                        pBNode->addValue(dof, mMGLevel, nodalVal);//外力に加算
                    }
                    break;
                case(ElementType::Prism2):
                    for(ivert=0; ivert < 15; ivert++){
                        integVal= pShPrism->getIntegValue15(ivert);
                        nodalVal= entVal * integVal;//等価節点力

                        pBNode= pBVol->getBNode(ivert);
                        pBNode->addValue(dof, mMGLevel, nodalVal);//外力に加算
                    }
                    break;
                default:
                    pLogger->Info(Utility::LoggerMode::Error, "CBoundaryVolumeMesh::distNeumannValue, invalid ElementType");
                    break;
            }//switch end
        };// idof ループ
    };//ivol ループ
}

//// Coarse Grid の最初の処理
////
//// Level==0:BNode周囲の要素集合平均値
////
//void CBoundaryVolumeMesh::distDirichletValue_at_CGrid()
//{
//    // 要素集合の平均をとるのはLevel=0の場合だけである
//    // --
//    if(mMGLevel!=0){
//        Utility::CLogger *pLogger= Utility::CLogger::Instance();
//        pLogger->Info(Utility::LoggerMode::Error, "MG-Level Error,  CBoundaryVolumeMesh::distDirichletValue");
//    }
//
//    //Solidの節点集合ごとに,DOFごとに,平均値をとる(Level==0)
//    //
//    uiint inode, numOfBNode=mvBNode.size();
//    CBoundaryNode *pBNode;
//    vuint vAggID; uiint numOfAggVol;
//
//    for(inode=0; inode < numOfBNode; inode++){
//        pBNode = mvBNode[inode];
//        vAggID = mvAggregateVol[inode];
//
//        numOfAggVol = vAggID.size();
//        uiint iVol, nVolID, nVolIndex;
//        CBoundaryVolume *pBVol;
//        double dCubicVol, dVal;//節点周囲のVolの体積、境界値
//        double dDirichletVal;  //自由度別の節点境界値
//        uiint idof, dof, numOfDOF=getNumOfDOF();
//
//        // 自由度別にディレクレ境界値を計算
//        //
//        for(idof=0; idof < numOfDOF; idof++){
//            dCubicVol=0.0; dVal=0.0;
//
//            //dof=pBNode->getDOF(idof);//自由度番号 <<<< VolumeMeshのBNodeの自由度番号管理とVolの自由度管理一致させる
//            dof=getDOF(idof);
//
//            for(iVol=0; iVol < numOfAggVol; iVol++){
//                nVolID = vAggID[iVol];
//                nVolIndex = mmBVolumeID2Index[nVolID];
//
//                pBVol = mvBVolume[nVolIndex];
//
//                dVal  += pBVol->getCubicVolume() * pBVol->getBndValue(dof);
//                dCubicVol += pBVol->getCubicVolume();
//            };
//            dDirichletVal = dVal/dCubicVol;// 節点値 = Σ(体積*値/Σ(体積)
//
//            pBNode->setValue(dof, mMGLevel, dDirichletVal);//境界値を代入
//        };
//    };
//}

//// Fine Grid
//// Level>=1:BNode間の平均
////
//void CBoundaryVolumeMesh::distDirichletValue_at_FGrid()
void CBoundaryVolumeMesh::distDirichletValue()
{
    //上位のグリッドの値をきめる
    //
    uiint ivol, numOfVol=mvBVolume.size();
    CBoundaryVolume *pBVol;
    for(ivol=0; ivol < numOfVol; ivol++){
        pBVol = mvBVolume[ivol];

        uiint idof, dof;
        for(idof=0; idof < getNumOfDOF(); idof++){
            dof = getDOF(idof);
            pBVol->distDirichletVal(dof, mMGLevel, mMaxMGLevel);
        };
    };
}


//// ディレクレ型境界値のBNodeへの再配分
////
//void CBoundaryVolumeMesh::distDirichletValue()
//{
//    Utility::CLogger *pLogger= Utility::CLogger::Instance();
//
//    // ディレクレ条件が設定されていなければ,Errorを返す.
//    // --
//    if(mnBndType != BoundaryType::Dirichlet){
//        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryVolumeMesh::distDirichletValue");
//        return;
//    }
//
//    if(mMGLevel==0){
//        distDirichletValue_at_CGrid();
//        if(mMaxMGLevel > 0) distDirichletValue_at_FGrid();
//    }else{
//        distDirichletValue_at_FGrid();
//    }
//
//}


// Refine 後処理 : 辺-面 BNode vectorの解放
//
void CBoundaryVolumeMesh::deleteProgData()
{
    uiint ivol, nNumOfBVol = mvBVolume.size();

    for(ivol=0; ivol < nNumOfBVol; ivol++){
        mvBVolume[ivol]->deleteProgData();
    };

    vector<CBoundaryNode*>().swap(mvBEdgeBNode);
    vector<CBoundaryNode*>().swap(mvBFaceBNode);
    vector<CBoundaryNode*>().swap(mvBVolBNode);
}







