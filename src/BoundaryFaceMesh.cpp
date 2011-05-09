//
//  BoundaryFaceMesh.cpp
//
//
//
//              2010.04.12
//              k.Takeda
#include <vector>

#include "BoundaryNode.h"
#include "BoundaryParts.h"
#include "BoundaryFace.h"

#include "BoundaryFaceMesh.h"
using namespace pmw;


CBoundaryFaceMesh::CBoundaryFaceMesh()
{
    ;
}

CBoundaryFaceMesh::~CBoundaryFaceMesh()
{
    for_each(mvBFace.begin(), mvBFace.end(), DeleteObject());
}


void CBoundaryFaceMesh::setBFace(const uint& index, CBoundaryFace* pBFace)
{
    uint id;
    id= pBFace->getID();
    
    mmBFaceID2Index[id]= index;

    mvBFace[index]= pBFace;
}
void CBoundaryFaceMesh::addBFace(CBoundaryFace* pBFace)
{
    mvBFace.push_back(pBFace);

    uint id;
    id= pBFace->getID();

    mmBFaceID2Index[id]= mvBFace.size()-1;
}


// 面-集合
// ----
void CBoundaryFaceMesh::resizeAggFace()
{
    uint numOfAgg= mvBNode.size();
    mvAggregateFace.resize(numOfAgg);// BNodeの面-集合

    //cout << "BoundaryFaceMesh::resizeAggFace, numOfAgg " << numOfAgg << endl;
}
// ----
// 点-周囲の面集合
// ----
void CBoundaryFaceMesh::setupAggFace()
{
    int numOfBNode= mvBNode.size();
    CBoundaryNode *pBNode;
    uint ibnode;
    //下位で収集した要素IDをクリア
    //----
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= mvBNode[ibnode];

        pBNode->clearAggElemID();
        pBNode->clearNeibElemVert();
    };
    
    uint numOfFace= mvBFace.size();
    CBoundaryFace *pBFace;
    uint iface;
    //BNodeに要素IDを集める
    //----
    for(iface=0; iface < numOfFace; iface++){
        pBFace= mvBFace[iface];
        
        //uint numOfVert= pBFace->getNumOfBNode();
        uint numOfVert= pBFace->getNumOfVert();
        uint ivert, faceID;
        for(ivert=0; ivert < numOfVert; ivert++){
            pBNode= pBFace->getBNode(ivert);
            
            faceID= pBFace->getID();
            pBNode->setAggElemID(faceID);
        };
    };

    //メンバー:AggreageFaceに移し替え
    //----
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= mvBNode[ibnode];

        uint numOfAgg= pBNode->getNumOfAggElem();
        uint iagg;

        for(iagg=0; iagg < numOfAgg; iagg++){
            mvAggregateFace[ibnode].push_back(pBNode->getAggElemID(iagg));
        };
    };
}
// ----
// 辺-周囲の面集合 => EdgeBNodeの生成
// ----
void CBoundaryFaceMesh::GeneEdgeBNode()
{
    uint countID= mvBNode.size();
    
    uint numOfBFace= mvBFace.size();
    uint iface;
    CBoundaryFace *pBFace;
    bool bfind(false);
    for(iface=0; iface < numOfBFace; iface++){
        pBFace= mvBFace[iface];
        
        uint numOfEdge= pBFace->getNumOfEdge();
        PairBNode pairBNode;
        uint iedge;
        
        for(iedge=0; iedge < numOfEdge; iedge++){
            if(!pBFace->isMarkingEdge(iedge)){
                pairBNode= pBFace->getPairBNode(iedge);

                uint numOfiAgg= pairBNode.first->getNumOfAggElem();
                uint numOfjAgg= pairBNode.second->getNumOfAggElem();
                uint iAgg,jAgg;
                uint iFaceID,jFaceID;
                bfind= false;

                // 隣接Faceがあるとして処理 => 無い場合は"bfind(偽)"
                // ----
                for(iAgg=0; iAgg < numOfiAgg; iAgg++){
                    iFaceID= pairBNode.first->getAggElemID(iAgg);

                    for(jAgg=0; jAgg < numOfjAgg; jAgg++){
                        jFaceID= pairBNode.second->getAggElemID(jAgg);

                        // 両端のBNodeの所有FaceIDが一致
                        // --------
                        if(iFaceID == jFaceID){
                            if(iFaceID != pBFace->getID()){

                                bfind= true;
                                // 生成
                                CBoundaryNode *pEdgeBNode= new CBoundaryNode;
                                
                                pEdgeBNode->setID(countID);
                                countID++;
                                
                                // FaceMeshにEdgeBNodeセット
                                mvBEdgeBNode.push_back(pEdgeBNode);

                                // pBFaceに"EdgeBNode,NeibFaceID"をセット
                                pBFace->setEdgeNeibFace(iedge,iFaceID);
                                pBFace->markingEdge(iedge);
                                pBFace->setEdgeBNode(iedge, pEdgeBNode);

                                // NeibFaceを取得
                                uint neib_index= mmBFaceID2Index[iFaceID];
                                CBoundaryFace* pNeibBFace= mvBFace[neib_index];
                                uint jedge= pNeibBFace->getEdgeID(pairBNode);

                                // NeibFaceに"EdgeBNode,pBFace-ID"をセット
                                pNeibBFace->setEdgeNeibFace(jedge, pBFace->getID());
                                pNeibBFace->markingEdge(jedge);
                                pNeibBFace->setEdgeBNode(jedge, pEdgeBNode);

                                //2次要素の場合 => 自身のmvBNodeにEdgeBNodeを追加
                                //----
                                if(pBFace->getOrder()==ElementOrder::Second){
                                    mvBNode.push_back(pEdgeBNode);
                                    mmBNodeID2Index[pEdgeBNode->getID()] = mvBNode.size()-1;
                                    
                                    //辺BNodeのAggFace
                                    uint index = mvBNode.size()-1;
                                    mvAggregateFace.resize(mvBNode.size());
                                    mvAggregateFace[index].push_back(pBFace->getID());
                                    mvAggregateFace[index].push_back(pNeibBFace->getID());

                                    //cout << "BoundaryFaceMesh::GeneEdgeBNode, 隣接Faceあり , mMGLevel=" << mMGLevel << endl;

                                    pBFace->replaceEdgeBNode();    //2次要素の辺Nodeを要素内mvBNodeに移設
                                    pNeibBFace->replaceEdgeBNode();//2次要素の辺Nodeを要素内mvBNodeに移設

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
                        
                    };//jAgg loop
                    if(bfind) break;
                    
                };//iAgg loop

                
                // 隣接Faceが存在しない場合
                // ----
                if(!bfind){
                    // 生成
                    CBoundaryNode *pEdgeBNode= new CBoundaryNode;
                    
                    pEdgeBNode->setID(countID);
                    countID++;
                    
                    // FaceMeshにEdgeBNodeセット
                    mvBEdgeBNode.push_back(pEdgeBNode);
                    
                    // pBFaceに"EdgeBNode"をセット
                    pBFace->setEdgeBNode(iedge, pEdgeBNode);
                    pBFace->markingEdge(iedge);

                    //2次要素の場合 => 自身のmvBNodeにEdgeBNodeを追加
                    //----
                    if(pBFace->getOrder()==ElementOrder::Second){
                        mvBNode.push_back(pEdgeBNode);
                        mmBNodeID2Index[pEdgeBNode->getID()] = mvBNode.size()-1;
                        
                        //辺BNodeのAggFace
                        uint index = mvBNode.size()-1;
                        mvAggregateFace.resize(mvBNode.size());
                        mvAggregateFace[index].push_back(pBFace->getID());

                        //cout << "BoundaryFaceMesh::GeneEdgeBNode, Non-Neib Face , mMGLevel=" << mMGLevel << endl;
                        pBFace->replaceEdgeBNode();//2次要素の辺Nodeを要素内mvBNodeに移設

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
                
            }//if(isMarkingEdge)
        };//iedge loop
        
        //BNodeにNodeをセット(全ての辺)
        pBFace->setupNode_Edge();

    };//iface loop
}


// FaceBNodeの生成
// ----
void CBoundaryFaceMesh::GeneFaceBNode()
{
    uint numOfFace= mvBFace.size();
    mvBFaceBNode.reserve(numOfFace);
    
    CBoundaryFace *pBFace;
    uint iface;
    uint countID= mvBNode.size() + mvBEdgeBNode.size();
    for(iface=0; iface < numOfFace; iface++){
        
        CBoundaryNode *pFaceBNode= new CBoundaryNode;
        
        // FaceMeshへFaceBNodeをセット
        mvBFaceBNode.push_back(pFaceBNode);
        
        // FaceBNodeの属性設定
        pFaceBNode->setMGLevel(mMGLevel+1);
        pFaceBNode->resizeValue(mMaxMGLevel-mMGLevel);
        pFaceBNode->setID(countID);
        countID++;
        
        // FaceにFaceBNodeをセット
        pBFace= mvBFace[iface];
        pBFace->setFaceBNode(pFaceBNode);
        
        //BNodeにNodeをセット(面)
        pBFace->setupNode_Face();
    };
}



// 1. 境界要素の再分割 => progBMeshにセット
// 2. BNode,EdgBNode,FaceBNode => progBMeshにセット
// ----
void CBoundaryFaceMesh::refine(CBoundaryFaceMesh *pProgBFaceMesh)
{
    CBoundaryFace *pBFace;
    uint iface, numOfFace= mvBFace.size();
    uint countID(0);//新ID
    vector<CBoundaryFace*> vprogBFace;//pBFaceから生成される子BFace
    //----
    //BFaceの分割 => progFaceMeshにセット
    //----
    for(iface=0; iface < numOfFace; iface++){
        pBFace= mvBFace[iface];

        pBFace->refine(countID, mvDOF);//Refine, 内部でcountID++
        
        vprogBFace.clear();
        vprogBFace= pBFace->getProgParts();
        
        uint i;
        for(i=0; i < vprogBFace.size(); i++){
            pProgBFaceMesh->addBFace(vprogBFace[i]);
        };
    };

    //----
    //BNode(頂点,辺,面) => progFaceMeshにセット
    //----
    uint numOfBNode    = mvBNode.size();
    //uint numOfEdgeBNode= mvBEdgeBNode.size();
    uint numOfFaceBNode= mvBFaceBNode.size();
    uint numOfProgBNode= numOfBNode + mnEdgeNodeCount + numOfFaceBNode;
    
    pProgBFaceMesh->resizeBNode(numOfProgBNode);
    
    uint ibnode;
    
    uint init= 0;
    uint end = numOfBNode;

    for(ibnode=init; ibnode < end; ibnode++){
        pProgBFaceMesh->setBNode(ibnode, mvBNode[ibnode]);
    };

    init= end;
    end = init + mnEdgeNodeCount;

    for(ibnode=init; ibnode < end; ibnode++){
        pProgBFaceMesh->setBNode(ibnode, mvBEdgeBNode[ibnode - init]);
    };

    init= end;
    end = init + numOfFaceBNode;

    for(ibnode=init; ibnode < end; ibnode++){
        pProgBFaceMesh->setBNode(ibnode, mvBFaceBNode[ibnode - init]);
    }
}




// ----
// Neumann条件の節点分配(形状関数による等価節点力:EquivalentNodeForce)
// ----
void CBoundaryFaceMesh::distNeumannValue()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    // ノイマン条件が設定されていなければ,Errorを返す.
    // --
    if(mnBndType != BoundaryType::Neumann){
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryFaceMesh::distNeumannValue");
        return;
    }
    
    // 形状関数
    CShapeQuad     *pQuadShape= CShapeQuad::Instance();
    CShapeTriangle *pTriShape= CShapeTriangle::Instance();

    CBoundaryNode *pBndNode;
    CBoundaryFace *pBndFace;
    double entVal;  //面の所有する境界値
    double integVal;//形状関数の積分値(正規化)
    double nodalVal;
    uint ivert, inode, numOfBNode= mvBNode.size();
    uint iface, numOfFace= mvBFace.size();
    uint idof, numOfDOF, dof;

    //節点力-ゼロクリア
    for(inode=0; inode < numOfBNode; inode++){
        pBndNode= mvBNode[inode];
        
        //numOfDOF= pBndNode->getNumOfDOF();
        numOfDOF = getNumOfDOF();
        for(idof=0; idof < numOfDOF; idof++){
            //dof= pBndNode->getDOF(idof);
            dof = getDOF(idof);
            pBndNode->initValue(dof, mMGLevel);//mgLevel別-自由度別の値をゼロクリア
        };
    };

    //節点力の計算-節点へ加算
    for(iface=0; iface < numOfFace; iface++){
        pBndFace= mvBFace[iface];

        //numOfDOF= pBndFace->getNumOfDOF();
        numOfDOF = getNumOfDOF();
        
        for(idof=0; idof < numOfDOF; idof++){

            dof = getDOF(idof);//自由度番号
            entVal= pBndFace->getBndValue(dof);//面の所有する値

            ////debug
            //cout << "CBoundaryFaceMesh::distNeumannValue, mgLevel=" << mMGLevel << ", dof=" << dof << ", entVal=" << entVal << endl;

            switch(pBndFace->getBFaceShape()){
                case(ElementType::Quad):
                    for(ivert=0; ivert < 4; ivert++){
                        // 等価節点力計算
                        integVal= pQuadShape->getIntegValue4(ivert);
                        nodalVal= entVal * integVal;

                        ////debug
                        //cout << "CBoundaryFaceMesh::distNeumannValue, mgLevel=" << mMGLevel << ", dof=" << dof << ", nodalVal=" << nodalVal << endl;

                        // 境界Node
                        pBndNode= pBndFace->getBNode(ivert);

                        // 境界値として節点"外力に加算"
                        pBndNode->addValue(dof, mMGLevel, nodalVal);
                    };
                    break;
                case(ElementType::Quad2):
                    for(ivert=0; ivert < 8; ivert++){
                        // 等価節点力計算
                        integVal= pQuadShape->getIntegValue8(ivert);
                        nodalVal= entVal * integVal;

                        // 境界Node
                        pBndNode= pBndFace->getBNode(ivert);

                        // 境界値として節点"外力に加算"
                        pBndNode->addValue(dof, mMGLevel, nodalVal);
                    };
                    break;
                case(ElementType::Triangle):
                    for(ivert=0; ivert < 3; ivert++){
                        // 等価節点力計算
                        integVal= pTriShape->getIntegValue3(ivert);
                        nodalVal= entVal * integVal;

                        // 境界Node
                        pBndNode= pBndFace->getBNode(ivert);

                        // 境界値として節点"外力に加算"
                        pBndNode->addValue(dof, mMGLevel, nodalVal);
                    };
                    break;
                case(ElementType::Triangle2):
                    for(ivert=0; ivert < 6; ivert++){
                        // 等価節点力計算
                        integVal= pTriShape->getIntegValue6(ivert);
                        nodalVal= entVal * integVal;

                        // 境界Node
                        pBndNode= pBndFace->getBNode(ivert);

                        // 境界値として節点"外力に加算"
                        pBndNode->addValue(dof, mMGLevel, nodalVal);
                    };
                    break;
                default:
                    //TODO:pLogger
                    break;
            }
        };
    };
}


// ----
// Dirichlet境界値のBNodeへの配分(Level==0の場合の最初の処理)
// ----
void CBoundaryFaceMesh::distDirichletValue_at_CGrid()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    // ディレクレ条件が設定されていなければ,Errorを返す.
    // --
    if(mnBndType != BoundaryType::Dirichlet){
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryFaceMesh::distDirichletValue");
        return;
    }
    
    //節点の面集合ごとに平均値をとる
    uint numOfBNode = mvBNode.size();
    uint inode;
    for(inode=0; inode < numOfBNode; inode++){

        CBoundaryNode *pBNode = mvBNode[inode];

        ////debug
        //cout << "CBoundaryFaceMesh::distDirichletValue_at_CGrid, mgLevel=" << mMGLevel
        //        << ", mvAggregateFace数=" << mvAggregateFace.size() << endl;

        vuint vAggFace = mvAggregateFace[inode];
        
        uint numOfAggFace = vAggFace.size();
        uint iface, nFaceID, nFaceIndex;
        CBoundaryFace *pBFace;
        double dArea, dVal, dDirichletVal;
        uint idof, numOfDOF=getNumOfDOF();
        uint dof;//自由度番号(idofは自由度配列のインデックス番号)


        // 自由度別にディレクレ境界値を計算
        //
        for(idof=0; idof < numOfDOF; idof++){

            dArea=0.0; dVal=0.0;
            dof = getDOF(idof);

            for(iface=0; iface < numOfAggFace; iface++){
                nFaceID = vAggFace[iface];
                nFaceIndex = mmBFaceID2Index[nFaceID];

                //cout << "CBoundaryFaceMesh::distDirichletValue_at_CGrid, FaceID " << nFaceID << endl;

                pBFace = mvBFace[nFaceIndex];

                dVal  += pBFace->getArea() * pBFace->getBndValue(dof);
                dArea += pBFace->getArea();
            };
            ////debug
            //cout << "CBoundaryFaceMesh::distDirichletValue_at_CGrid, mgLevel=" << mMGLevel
            //        << ", dVal=" << dVal << ", dArea=" << dArea << endl;

            dDirichletVal = dVal/dArea;// 節点値 = Σ(面積*値)/Σ(面積)

            pBNode->setValue(dof, mMGLevel, dDirichletVal);//境界値を代入
        };
    };
}

// Fine Grid
// Level>=1:BNode間の平均
//
void CBoundaryFaceMesh::distDirichletValue_at_FGrid()
{
    //上位のグリッドの値をきめる
    //
    uint iface, numOfFace=mvBFace.size();
    CBoundaryFace *pBFace;
    for(iface=0; iface < numOfFace; iface++){
        pBFace = mvBFace[iface];

        uint idof, dof;
        for(idof=0; idof < getNumOfDOF(); idof++){
            dof = getDOF(idof);
            pBFace->distDirichletVal(dof, mMGLevel, mMaxMGLevel);
        };
    };
}


// Refine時のデータの解放
// ----
// 各Faceの辺BNodeの解放
// mvBEdgeBNodeの解放
// mvBFaceBNodeの解放
void CBoundaryFaceMesh::deleteProgData()
{
    uint iface, nNumOfBFace=mvBFace.size();

    for(iface=0; iface < nNumOfBFace; iface++) mvBFace[iface]->deleteProgData();

    vector<CBoundaryNode*>().swap(mvBEdgeBNode);
    vector<CBoundaryNode*>().swap(mvBFaceBNode);
}






