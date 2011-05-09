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


void CBoundaryVolumeMesh::resizeVolume(const uint& res_size)
{
    mvBVolume.resize(res_size);
}


void CBoundaryVolumeMesh::setBVolume(const uint& index, CBoundaryVolume* pBVolume)
{
    uint id;
    id= pBVolume->getID();

    mmBVolumeID2Index[id]= index;

    mvBVolume[index]= pBVolume;
}


void CBoundaryVolumeMesh::addBVolume(CBoundaryVolume* pBVolume)
{
    mvBVolume.push_back(pBVolume);

    uint id;
    id= pBVolume->getID();

    mmBVolumeID2Index[id]= mvBVolume.size()-1;
}



// 頂点のVolID集合の領域確保
// ----
void CBoundaryVolumeMesh::resizeAggVol()
{
    uint numOfAgg= mvBNode.size();
    
    mvAggregateVol.resize(numOfAgg);
}

// 頂点へVolumeIDを集合
// ----
void CBoundaryVolumeMesh::setupAggVol()
{
    uint numOfBNode= mvBNode.size();
    CBoundaryNode  *pBNode;
    uint ibnode;
    //下位グリッドのAggデータをクリア
    //----
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= mvBNode[ibnode];

        pBNode->clearAggElemID();
        pBNode->clearNeibElemVert();
    };
    
    uint numOfVol= mvBVolume.size();
    CBoundaryVolume *pBVol;
    uint ivol;
    //BNodeの周囲Volを集める
    //----
    for(ivol=0; ivol < numOfVol; ivol++){
        pBVol= mvBVolume[ivol];
        
        uint numOfVert= pBVol->getNumOfBNode();
        uint ivert;
        for(ivert=0; ivert < numOfVert; ivert++){
            pBNode= pBVol->getBNode(ivert);
            pBNode->setAggElemID(pBVol->getID());
        };
    };
    
    //Aggregateに移す.
    //----
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= mvBNode[ibnode];

        uint numOfAgg= pBNode->getNumOfAggElem();
        uint iAgg;

        for(iAgg=0; iAgg < numOfAgg; iAgg++){
            mvAggregateVol[ibnode].push_back(pBNode->getAggElemID(iAgg));
        };
    };
}


// EdgeBNode生成
//
void CBoundaryVolumeMesh::GeneEdgeBNode()
{
    uint numOfVol= mvBVolume.size();
    uint ivol;
    CBoundaryVolume *pBVol;

    CElement *pElem;//各BVolumeが所有するElement
    CNode *pNode0,*pNode1;//辺の両端のNode(MeshのNode)
    CNode *pEdgeNode;//辺中央のNode

    uint countID= mvBNode.size();//EdgeBNodeのID初期値

    for(ivol=0; ivol < numOfVol; ivol++){
        pBVol= mvBVolume[ivol];

        uint numOfEdge= pBVol->getNumOfEdge();
        uint iedge;
        PairBNode pairBNode;
        for(iedge=0; iedge < numOfEdge; iedge++){

            if(!pBVol->isMarkingEdge(iedge)){

                bool bfind(false);//辺に隣接Volが存在するか否かの真偽

                pairBNode= pBVol->getPairBNode(iedge);

                uint numOfiAgg= pairBNode.first->getNumOfAggElem();
                uint numOfjAgg= pairBNode.second->getNumOfAggElem();
                uint iAgg,jAgg, iVolID, jVolID;
                for(iAgg=0; iAgg < numOfiAgg; iAgg++){
                    iVolID= pairBNode.first->getAggElemID(iAgg);

                    for(jAgg=0; jAgg < numOfjAgg; jAgg++){
                        jVolID= pairBNode.second->getAggElemID(jAgg);

                        if(iVolID==jVolID){
                            if(iVolID != pBVol->getID()){

                                bfind= true;

                                // 生成
                                CBoundaryNode *pEdgeBNode= new CBoundaryNode;// <<<<<<<<<<<<<< new

                                // EdgeBNodeの属性設定
                                pEdgeBNode->setMGLevel(mMGLevel+1);
                                pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel);
                                pEdgeBNode->setID(countID);
                                countID++;

                                // VolumeMeshにEdgeBNodeセット
                                mvBEdgeBNode.push_back(pEdgeBNode);

                                // pBVolの"辺"にEdgeBNodeをセット
                                pBVol->setEdgeNeibVol(iedge, iVolID);
                                pBVol->markingEdge(iedge);
                                pBVol->setEdgeBNode(iedge, pEdgeBNode);


                                // 隣接するVolumeを取得
                                uint neibIndex= mmBVolumeID2Index[iVolID];
                                CBoundaryVolume *pNeibBVol= mvBVolume[neibIndex];
                                uint jedge= pNeibBVol->getEdgeID(pairBNode);

                                // pNeibBVol(隣接Vol)の"辺"にEdgeBNodeをセット
                                pNeibBVol->setEdgeNeibVol(jedge, pBVol->getID());
                                pNeibBVol->markingEdge(jedge);
                                pNeibBVol->setEdgeBNode(jedge, pEdgeBNode);
                                

                                //辺中央のNodeを取得してEdgeBNodeにセット
                                //----
                                pElem = pBVol->getElement();
                                pNode0 = pairBNode.first->getNode();
                                pNode1 = pairBNode.second->getNode();
                                uint nEdgeIndex = pElem->getEdgeIndex(pNode0, pNode1);
                                pEdgeNode = pElem->getEdgeInterNode(nEdgeIndex);
                                pEdgeBNode->setNode(pEdgeNode);
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

                    // EdgeBNodeの属性設定
                    pEdgeBNode->setMGLevel(mMGLevel+1);
                    pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel);
                    pEdgeBNode->setID(countID);
                    countID++;

                    // VolumeMeshにEdgeBNodeセット
                    mvBEdgeBNode.push_back(pEdgeBNode);

                    // pBVolの"辺"にEdgeBNodeをセット
                    pBVol->markingEdge(iedge);
                    pBVol->setEdgeBNode(iedge, pEdgeBNode);


                    //辺中央のNodeを取得してEdgeBNodeにセット
                    //----
                    pElem = pBVol->getElement();
                    pNode0 = pairBNode.first->getNode();
                    pNode1 = pairBNode.second->getNode();
                    uint nEdgeIndex = pElem->getEdgeIndex(pNode0, pNode1);
                    pEdgeNode = pElem->getEdgeInterNode(nEdgeIndex);
                    pEdgeBNode->setNode(pEdgeNode);
                }
            }// if(isMarkingEdge)
        };// iedge loop
    };// ivol loop
}

// FaceBNode生成
//
void CBoundaryVolumeMesh::GeneFaceBNode()
{
    uint numOfVol= mvBVolume.size();
    CBoundaryVolume *pBVol;
    CBoundaryVolume *pNeibBVol;

    CElement *pElem;//各BVolumeが所有するElement
    vector<CNode*> vNode; vNode.resize(3);//面構成Node(MeshのNode)
    CNode *pFaceNode;//面中央のNode

    uint countID= mvBNode.size() + mvBEdgeBNode.size();//FaceBNodeのID初期値
    uint ivol;
    
    // 体積
    for(ivol=0; ivol < numOfVol; ivol++){
        
        pBVol= mvBVolume[ivol];
        
        uint numOfFace= pBVol->getNumOfFace();
        uint iface;
        vector<CBoundaryNode*> vBNode;

        // 面構成BNode
        for(iface=0; iface < numOfFace; iface++){

            if(!pBVol->isMarkingFace(iface)){

                vBNode= pBVol->getFaceCnvNodes(iface);

                uint iVolID, jVolID, kVolID;

                uint iAgg, jAgg, kAgg;
                uint numOfiAgg= vBNode[0]->getNumOfAggElem();
                uint numOfjAgg= vBNode[1]->getNumOfAggElem();
                uint numOfkAgg= vBNode[2]->getNumOfAggElem();

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
                                uint neibIndex= mmBVolumeID2Index[iVolID];
                                pNeibBVol= mvBVolume[neibIndex];
                                uint neibFace= pNeibBVol->getFaceID(vBNode);

                                pNeibBVol->setFaceNeibVol(neibFace, pBVol->getID());
                                pNeibBVol->markingFace(neibFace);
                                pNeibBVol->setFaceBNode(neibFace, pFaceBNode);

                                
                                // Nodeのセット
                                // Elementから面Nodeを取得して、FaceBNodeにセット
                                // ----
                                pElem = pBVol->getElement();
                                uint ivert;
                                for(ivert=0; ivert < 3; ivert++){
                                    vNode[ivert] = vBNode[ivert]->getNode();
                                };
                                uint nFaceIndex = pElem->getFaceIndex(vNode[0],vNode[1],vNode[2]);
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
                    uint ivert;
                    for(ivert=0; ivert < 3; ivert++){
                        vNode[ivert] = vBNode[ivert]->getNode();
                    };
                    uint nFaceIndex = pElem->getFaceIndex(vNode[0],vNode[1],vNode[2]);
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
    uint numOfVol= mvBVolume.size();
    CBoundaryVolume *pBVol;
    CElement *pElem;//各BVoluemが所有するElement
    CNode *pNode;//Element中央のNode

    uint countID= mvBNode.size() + mvBEdgeBNode.size() + mvBFaceBNode.size();//VolumeBNodeのID初期値
    uint ivol;

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
    uint ivol, numOfVol= mvBVolume.size();
    vector<CBoundaryVolume*> vProgVol;// pBVolから生成される子BVolume
    //----
    //BVolumeの分割 => progBVolMeshにセット
    //----
    uint countID(0);//新ID
    for(ivol=0; ivol < numOfVol; ivol++){
        pBVol= mvBVolume[ivol];

        pBVol->refine(countID, mvDOF);//  Refine, 内部でcountID++

        vProgVol.clear();
        vProgVol= pBVol->getProgParts();

        uint i;
        for(i=0; i < vProgVol.size(); i++){
            pProgBVolMesh->addBVolume(vProgVol[i]);
        };
    };

    //----
    //BNode(頂点,辺,面) => progBVolMeshにセット
    //----
    uint numOfBNode    = mvBNode.size();
    uint numOfEdgeBNode= mvBEdgeBNode.size();
    uint numOfFaceBNode= mvBFaceBNode.size();
    uint numOfVolBNode = mvBVolBNode.size();
    uint numOfProgBNode= numOfBNode + numOfEdgeBNode + numOfFaceBNode + numOfVolBNode;

    pProgBVolMesh->resizeBNode(numOfProgBNode);
    
    uint ibnode;
    // 頂点BNode
    uint init = 0;
    uint end  = numOfBNode;
    
    for(ibnode=init; ibnode < end; ibnode++){
        pProgBVolMesh->setBNode(ibnode, mvBNode[ibnode]);
    };
    
    // 辺BNode
    init = numOfBNode;
    end  = numOfBNode + numOfEdgeBNode;
    
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

    uint inode, numOfBNode=mvBNode.size();
    CBoundaryNode *pBNode;
    uint dof, idof;
    //節点力-ゼロクリア
    for(inode=0; inode < numOfBNode; inode++){
        pBNode= mvBNode[inode];

        //uint numOfDOF= pBNode->getNumOfDOF();
        uint numOfDOF = getNumOfDOF();

        for(idof=0; idof < numOfDOF; idof++){
            //dof= pBNode->getDOF(idof);
            dof = getDOF(idof);
            pBNode->initValue(dof, mMGLevel);//mgLevel別-自由度別の値をゼロクリア
        };
    };

    //debug
    cout << "CBoundaryVolumeMesh::distNeumannValue, mgLevel=" << mMGLevel << ", aaaaaaaaaaa" << endl;

    //形状関数の積分値による配分
    double integVal, nodalVal, entVal;
    uint ivol, numOfVol=mvBVolume.size();
    uint numOfDOF;
    for(ivol=0; ivol < numOfVol; ivol++){
        CBoundaryVolume *pBVol = mvBVolume[ivol];

        numOfDOF = getNumOfDOF();

        //debug
        cout << "CBoundaryVolumeMesh::distNeumannValue, mgLevel=" << mMGLevel << ", numOfDOF=" << numOfDOF << endl;

        // 自由度別の外力配分
        // ----
        for(idof=0; idof < numOfDOF; idof++){

            dof = getDOF(idof);//自由度番号
            entVal = pBVol->getBndValue(dof);//体積の所有する境界値

            //debug
            cout << "CBoundaryVolumeMesh::distNeumannValue, mgLevel=" << mMGLevel << ", dof=" << dof << ", entVal=" << entVal << endl;

            uint ivert;
            switch(pBVol->getElemType()){
                case(ElementType::Hexa):
                    for(ivert=0; ivert < 8; ivert++){
                        integVal= pShHexa->getIntegralValue8(ivert);
                        nodalVal= entVal * integVal;// 等価節点力計算
                        
                        //debug
                        cout << "CBoundaryVolumeMesh::distNeumannValue, mgLevel=" << mMGLevel << ", dof=" << dof << ", nodalVal=" << nodalVal << endl;
                        

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

// Coarse Grid の最初の処理
//
// Level==0:BNode周囲の要素集合平均値
//
void CBoundaryVolumeMesh::distDirichletValue_at_CGrid()
{
    // 要素集合の平均をとるのはLevel=0の場合だけである
    // --
    if(mMGLevel!=0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "MG-Level Error,  CBoundaryVolumeMesh::distDirichletValue");
    }

    //Solidの節点集合ごとに,DOFごとに,平均値をとる(Level==0)
    //
    uint inode, numOfBNode=mvBNode.size();
    CBoundaryNode *pBNode;
    vuint vAggID; uint numOfAggVol;

    for(inode=0; inode < numOfBNode; inode++){
        pBNode = mvBNode[inode];
        vAggID = mvAggregateVol[inode];

        numOfAggVol = vAggID.size();
        uint iVol, nVolID, nVolIndex;
        CBoundaryVolume *pBVol;
        double dCubicVol, dVal;//節点周囲のVolの体積、境界値
        double dDirichletVal;  //自由度別の節点境界値
        uint idof, dof, numOfDOF=getNumOfDOF();

        // 自由度別にディレクレ境界値を計算
        //
        for(idof=0; idof < numOfDOF; idof++){
            dCubicVol=0.0; dVal=0.0;

            //dof=pBNode->getDOF(idof);//自由度番号 <<<< VolumeMeshのBNodeの自由度番号管理とVolの自由度管理一致させる
            dof=getDOF(idof);

            for(iVol=0; iVol < numOfAggVol; iVol++){
                nVolID = vAggID[iVol];
                nVolIndex = mmBVolumeID2Index[nVolID];

                pBVol = mvBVolume[nVolIndex];

                dVal  += pBVol->getCubicVolume() * pBVol->getBndValue(dof);
                dCubicVol += pBVol->getCubicVolume();
            };
            dDirichletVal = dVal/dCubicVol;// 節点値 = Σ(体積*値/Σ(体積)

            pBNode->setValue(dof, mMGLevel, dDirichletVal);//境界値を代入
        };
    };
}

// Fine Grid
// Level>=1:BNode間の平均
//
void CBoundaryVolumeMesh::distDirichletValue_at_FGrid()
{
    //上位のグリッドの値をきめる
    //
    uint ivol, numOfVol=mvBVolume.size();
    CBoundaryVolume *pBVol;
    for(ivol=0; ivol < numOfVol; ivol++){
        pBVol = mvBVolume[ivol];

        uint idof, dof;
        for(idof=0; idof < getNumOfDOF(); idof++){
            dof = getDOF(idof);
            pBVol->distDirichletVal(dof, mMGLevel);
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
    uint ivol, nNumOfBVol = mvBVolume.size();

    for(ivol=0; ivol < nNumOfBVol; ivol++){
        mvBVolume[ivol]->deleteProgData();
    };

    vector<CBoundaryNode*>().swap(mvBEdgeBNode);
    vector<CBoundaryNode*>().swap(mvBFaceBNode);
    vector<CBoundaryNode*>().swap(mvBVolBNode);
}







