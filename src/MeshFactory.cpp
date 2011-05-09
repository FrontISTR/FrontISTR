//
//  MeshFactory.cpp
//
//
//
//                      2009.09.29
//			2008.11.05
//			k.Takeda
#include "ElementType.h"
#include "BoundaryMesh.h"
#include "BoundaryParts.h"
#include "BoundaryFace.h"
#include "BndVertex.h"
#include "BoundaryNodeMesh.h"
#include "BoundaryFaceMesh.h"
#include "CommNode.h"
#include "CommMesh2.h"
#include "BoundaryNode.h"


#include <vector>

#include "FaceTree.h"
#include "EdgeTree.h"
#include "ContactNode.h"
#include "GMGModel.h"
#include "CommFace.h"
#include "ContactMesh.h"
#include "AssyModel.h"


#include "Element.h"


#include "Mesh.h"
#include "MeshFactory.h"
#include "BoundaryHexa.h"
#include "SkinFace.h"
using namespace pmw;

// construct & destruct
//
CMeshFactory::CMeshFactory(void)
{
    mpLogger = Utility::CLogger::Instance();
}

CMeshFactory::~CMeshFactory(void)
{
    cout << "~CMeshFactory" << endl;
}

// MGを構築しない場合は、MeshのsetupAggregateだけを構築
//
void CMeshFactory::refineMesh()
{
    if(mMGLevel > 0){
        MGMeshConstruct();// Multi Gridを構築
    }else{
        SGMeshConstruct();// Single Gridを構築
    }
}

// Single Grid
//
void CMeshFactory::SGMeshConstruct()
{
    CAssyModel *pAssy;
    CMesh      *pMesh;

    pAssy = mpGMGModel->getAssyModel(0);

    uint nLevel=0;
    uint numOfMesh = pAssy->getNumOfMesh();
    uint imesh;
    for(imesh=0; imesh < numOfMesh; imesh++){
        pMesh = pAssy->getMesh(imesh);

        pMesh->setupAggregate(nLevel);

        string str = "SGMeshConstruct , pMesh->setupAggElement  at " + boost::lexical_cast<string>(imesh);
        mpLogger->Info(Utility::LoggerMode::MWDebug, str);
    };
}

// Refine for MultiGrid
// --
// current Mesh      => pMesh
// prolongation Mesh => pProgMesh
// 
// current CommMesh  => pCommMesh (通信領域メッシュ)
// prolongation CommMesh => pProgCommMesh
// 
void CMeshFactory::MGMeshConstruct()
{
    // AseeyModel & Mesh --
    CAssyModel *pAssy, *pProgAssy;
    CMesh      *pMesh, *pProgMesh;
    // Element
    CElement         *pElem=NULL;//元になった要素
    vector<CElement*> vProgElem; //分割された新しい要素達

    // 通信Mesh
    CCommMesh  *pCommMesh, *pProgCommMesh;
    // 通信要素(CommElem)
    CCommElement      *pCommElem;//元になったCommElem(親のCommElem)
    vector<CCommElement*> vProgCommElem;//生成されるprogCommElemのコンテナ
    
    uint numOfCommMesh,numOfCommElemAll,numOfCommNode;
    uint icommesh,icomelem,iprocom;
    

    // 階層”0”のMesh数を基準に各階層のAssyModelのMeshを決める.
    // ---
    pAssy= mpGMGModel->getAssyModel(0);
    uint numOfMesh= pAssy->getNumOfMesh();
    
    uint ilevel,imesh,ielem;
    // ---
    // 階層Level ループ
    // ---
    for(ilevel=0; ilevel< mMGLevel; ilevel++){
        
        pAssy= mpGMGModel->getAssyModel(ilevel);
        
        // prolongation AssyModel
        //
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);//FileReadRefineブロックでAssyModelは生成済み
        pProgAssy->resizeMesh(numOfMesh);

        pProgAssy->intializeBucket(pAssy->getMaxMeshID(),pAssy->getMinMeshID());//逆引き配列の領域確保
        pProgAssy->setMaxMeshID(pAssy->getMaxMeshID());
        pProgAssy->setMinMeshID(pAssy->getMinMeshID());
        // ---
        // Mesh(パーツ) ループ in AssyMode
        // ---
        for(imesh=0; imesh< numOfMesh; imesh++){

            pMesh= pAssy->getMesh(imesh);//Current_Level Mesh(最初期はLevel==0：ファイル読み込み時)
            pProgAssy->setBucket(pMesh->getMeshID(), imesh);//progAssyの逆引きに"id-index"をセット
            
            // <<<< start ::pProgMeshの生成処理 >>>>
            //
            pProgMesh = new CMesh;          //Upper_Level Mesh ==(prolongation Mesh)
            pProgMesh->setMGLevel(ilevel+1);//上位MeshのMultiGridレベルを設定(初期pMeshはファイル読み込み時のLevel==0)
            pProgMesh->setMeshID(pMesh->getMeshID());//Mesh_ID は,同一のIDとする.

            
            // Refineの準備, 頂点集合の要素と,辺-面-体積中心の節点生成(progMeshの節点生成),辺-面の要素集合
            //
            if(ilevel==0){
                // 条件(ilevel >= 1) は, pProgMesh->setupAggregate()をMesh-Refine後に行ってあるので不要:CommMeshのRefineの為.
                pMesh->setupAggregate(ilevel); //Node集合Element, Node集合Nodeの計算
                mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupAggElement finish at ilevel==0");
            }
            pMesh->presetProgMesh(pProgMesh);//prolongation_Meshのノード,要素リザーブ(reserve) && pMeshのノード,要素をセット
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->presetProgMesh finish");
            pMesh->setupEdgeElement(pProgMesh, ilevel);//辺(Edge)節点, 辺に集合する要素の計算
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupEdgeElement finish");
            pMesh->setupFaceElement(pProgMesh);//面(Face)節点, 面に隣接する要素の計算
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupFaceElement finish");
            pMesh->setupVolumeNode(pProgMesh); //体(Volume)節点:要素中心の節点
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupVolumeNode finish");

            pMesh->replaceEdgeNode();//2次要素の場合、辺ノードを要素ノードとして移し替え
            
            // 新ElementID
            uint numOfElem = pMesh->getNumOfElement();
            uint elementID= 0;// 新たに生成される要素のIDを割振る. -> 各divid関数でカウントアップ.
                              // ElementのID(Index)初期値は,土台のMeshとは無関係 <= Nodeとは異なる.

            for(ielem=0; ielem< numOfElem; ielem++){
                pElem= pMesh->getElement(ielem);
                vProgElem.clear();
                
                GeneProgElem(ilevel, pElem, vProgElem, elementID, pProgMesh);//再分割要素の生成, 2010.5.31VC++同様に変更
            };
            
            // ノード, 数要素数のセット
            pProgMesh->setupNumOfNode();
            pProgMesh->setupNumOfElement();
            pProgMesh->setSolutionType(mnSolutionType);
            
            
            // pProgMeshの AggregateElement, AggregateNodeの生成
            // --
            uint numOfNode= pProgMesh->getNodeSize();
            CAggregateElement *pAggElem;
            CAggregateNode    *pAggNode;
            pProgMesh->resizeAggregate(numOfNode);

            if(mnSolutionType==SolutionType::FEM){
                for(uint iAgg=0; iAgg< numOfNode; iAgg++){
                    pAggElem = new CAggregateElement;
                    pProgMesh->setAggElement(pAggElem, iAgg);
                };
            }
            if(mnSolutionType==SolutionType::FVM){
                for(uint iAgg=0; iAgg< numOfNode; iAgg++){
                    pAggNode = new CAggregateNode;
                    pProgMesh->setAggNode(pAggNode, iAgg);
                };
            }
            // !注意 CMesh::setupAggregate()は,このルーチンの先頭のRefine準備でpMeshに対してコールするのでpProgMeshにはコールしない.
            
            // prolongation AssyModelに,pProgMeshをセット
            //
            pProgAssy->setMesh(pProgMesh,pMesh->getMeshID());//progAssyModel に,progMeshをセット
            

            // progMeshのNode逆引きセットアップ(CommMeshのRefineでprogMeshのBucketを使用)
            // --
            uint progNodeSize = pProgMesh->getNodeSize();
            pProgMesh->initBucketNode(progNodeSize+1, 0);//Bucketの領域確保
            pProgMesh->setupBucketNode();                //Bucketの"ID,Index"の一括処理
            // progMeshのElement逆引きセットアップ
            // --
            uint progElemSize = pProgMesh->getElementSize();
            pProgMesh->initBucketElement(progElemSize+1,0);
            pProgMesh->setupBucketElement();

            // <<<< end ::pProgMeshの生成処理 >>>>
            cout << "pProgMesh ノード数 == " << pProgMesh->getNumOfNode() << endl;//debug

            // Elementグループ
            //
            uint iGrp, nNumOfGrp = pMesh->getNumOfElemGrp();
            CElementGroup *pElemGrp, *pProgElemGrp;
            for(iGrp=0; iGrp < nNumOfGrp; iGrp++){

                pProgElemGrp = new CElementGroup;

                pElemGrp = pMesh->getElemGrpIX(iGrp);

                pElemGrp->refine(pProgElemGrp);

                pProgElemGrp->setMesh(pProgMesh);
                pProgElemGrp->setID(pElemGrp->getID());
                pProgElemGrp->setName(pElemGrp->getName());

                pProgMesh->addElemGrp(pProgElemGrp);
            };
            

            // progCommMeshの前処理
            // --
            pProgMesh->setupAggregate(ilevel+1);


            ////////////////////////////////////////////////
            // <<<< start ::pProgCommCommMeshの生成処理 >>>>
            //
            numOfCommMesh= pMesh->getNumOfCommMesh();
            // --
            // CommMesh(通信Mesh)ループ in Mesh
            // --
            for(icommesh=0; icommesh< numOfCommMesh; icommesh++){
                pCommMesh= pMesh->getCommMesh(icommesh);

                // "new CommMesh"に,下段階層CommMeshのプロパティをセット
                // --
                CIndexBucket *pBucket= pProgMesh->getBucket();
                pProgCommMesh= new CCommMesh(pBucket);// <<<<<<<<<<<<<<-- prolongation CommMesh
                pProgCommMesh->setCommID( pCommMesh->getCommID());
                pProgCommMesh->setRankID( pCommMesh->getRankID());
                pProgCommMesh->setTransmitRankID( pCommMesh->getTransmitRankID());

                pProgMesh->setCommMesh(pProgCommMesh);//progMeshにprogCommMeshをセット


                numOfCommElemAll= pCommMesh->getNumOfCommElementAll();//カレントCommElement数
                pProgCommMesh->reserveCommElementAll(numOfCommElemAll*8);// <<<-- CommElemAllリザーブ:下段階層のCommElem数の8倍

                numOfCommNode= pCommMesh->getNumOfNode();//カレントCommMeshのNode数
                pProgCommMesh->reserveNode(numOfCommNode*8);// <<<-- ノード・リザーブ:下段階層のCommMeshのノード数の8倍
                // --
                // CommElementループ in CommMesh (カレントCommMesh)
                // --
                for(icomelem=0; icomelem< numOfCommElemAll; icomelem++){
                    //pCommMeshからCommElemを取得し,progCommElemを生成 => progCommMeshにセット, rankはCommMeshが所有
                    //
                    pCommElem= pCommMesh->getCommElementAll(icomelem);
                    pCommElem->setupProgNodeRank(ilevel+1);//辺,面,体積中心にRank設定 <<<<-- progCommElemのNodeランク

                    vProgCommElem.clear();
                    GeneProgCommElem(pCommElem, vProgCommElem);//<<<<<<-- prolongation CommElementの生成

                    // --
                    // progCommMeshへprogCommElemAllのセット(全てのCommElement)
                    // --
                    for(iprocom=0; iprocom< vProgCommElem.size(); iprocom++){
                        pProgCommMesh->setCommElementAll(vProgCommElem[iprocom]);
                    };
                };//CommElementループ・エンド

                //1.CommMesh内でのCommElemのIndex番号の割り振り && CommMesh内のCommElementを通信するか否かのCommElementに選別
                pProgCommMesh->AllocateCommElement();
                
                //2.CommMesh内でのCommElemの隣接情報
                //3.CommMesh内でのNodeのIndex番号の割り振り,CommMeshのmvNode,mvSendNode,mvRecvNodeの取得
                pProgCommMesh->setupAggCommElement(pProgMesh->getElements());
                pProgCommMesh->sortCommNodeIndex();// CommMesh内でのNode Index番号生成, Send,Recvノードの選別, DNode,DElementの選別ソート
                                                   // mvNodeのセットアップもsortCommNodeIndexから,CommElementのsetCommNodeIndex内でmvNodeにセットしている.
                //4.mapデータのセットアップ
                pProgCommMesh->setupMapID2CommID();

            };//CommMeshループ・エンド
            

            // Mesh のNode,Elementの計算領域整理:MeshのmvNode,mvElementから計算に使用しないNode(DNode),Element(DElement)を移動
            // --
            pProgMesh->sortMesh();

            // Meshが,ソートされたので,Bucketを再セットアップ
            pProgMesh->setupBucketNode();//Bucketの"ID,Index"の一括処理
            pProgMesh->setupBucketElement();
            //
            // <<<< end ::pProgCommCommMeshの生成処理 >>>>

        };//imesh ループ エンド
    };//ilevel ループ エンド
    
    //
    // 2次要素の場合：最終LevelのMeshに辺ノードをprogMeshに生成
    //
    pAssy= mpGMGModel->getAssyModel(mMGLevel);//最終LevelのAssyModel
    
    numOfMesh= pAssy->getNumOfMesh();
    for(imesh=0; imesh< numOfMesh; imesh++){
        pMesh= pAssy->getMesh(imesh);
        
        pMesh->setupEdgeElement(NULL, mMGLevel);//2次要素として利用するため,最終レベルのMeshに辺ノードを生成
        pMesh->replaceEdgeNode();//辺ノードを移し替え
    };
}

// 再分割要素(progElem)の生成 2010.05.31VC++同様に変更
//
void CMeshFactory::GeneProgElem(const uint& ilevel,CElement* pElem, vector<CElement*>& vProgElem, uint& elementID, CMesh* pProgMesh)
{
    CElement *pProgElem;
    uint i;
    // divid Element(要素の分割)
    switch(pElem->getType()){
        case(ElementType::Hexa):
            vProgElem.reserve(8);//分割された新しい要素の生成
            for(i=0; i< 8; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();

                vProgElem.push_back(pProgElem);
            };
            dividHexa(pElem,vProgElem, elementID, pProgMesh);
            break;

        case(ElementType::Hexa2):
            vProgElem.reserve(8);//分割された新しい要素の生成
            for(i=0; i< 8; i++){
                pProgElem= new CHexa2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();

                vProgElem.push_back(pProgElem);
            };
            dividHexa(pElem,vProgElem, elementID, pProgMesh);
            break;

        case(ElementType::Tetra):
            vProgElem.reserve(4);//分割された新しい要素の生成
            for(i=0; i< 4; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();

                vProgElem.push_back(pProgElem);
            };
            dividTetra(pElem,vProgElem, elementID, pProgMesh);
            break;

        case(ElementType::Tetra2):
            vProgElem.reserve(4);//分割された新しい要素の生成
            for(i=0; i< 4; i++){
                pProgElem= new CHexa2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                
                vProgElem.push_back(pProgElem);
            };
            dividTetra(pElem,vProgElem, elementID, pProgMesh);
            break;

        case(ElementType::Prism):
            vProgElem.reserve(6);//分割された新しい要素の生成
            for(i=0; i< 6; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();

                vProgElem.push_back(pProgElem);
            };
            dividPrism(pElem,vProgElem, elementID,pProgMesh);
            break;

        case(ElementType::Prism2):
            vProgElem.reserve(6);//分割された新しい要素の生成
            for(i=0; i< 6; i++){
                pProgElem= new CHexa2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                
                vProgElem.push_back(pProgElem);
            };
            dividPrism(pElem,vProgElem, elementID,pProgMesh);
            break;

        case(ElementType::Quad):
            vProgElem.reserve(4);//分割された新しい要素の生成
            for(i=0; i< 4; i++){
                pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();

                vProgElem.push_back(pProgElem);
            };
            dividQuad(pElem,vProgElem, elementID,pProgMesh);
            break;

        case(ElementType::Quad2):
            vProgElem.reserve(4);//分割された新しい要素の生成
            for(i=0; i< 4; i++){
                pProgElem= new CQuad2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();

                vProgElem.push_back(pProgElem);
            };
            dividQuad(pElem,vProgElem, elementID,pProgMesh);
            break;

        case(ElementType::Triangle):
            vProgElem.reserve(3);//分割された新しい要素の生成
            for(i=0; i< 3; i++){
                pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();

                vProgElem.push_back(pProgElem);
            };
            dividTriangle(pElem,vProgElem, elementID,pProgMesh);
            break;

        case(ElementType::Triangle2):
            vProgElem.reserve(3);//分割された新しい要素の生成
            for(i=0; i< 3; i++){
                pProgElem= new CQuad2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();

                vProgElem.push_back(pProgElem);
            };
            dividTriangle(pElem,vProgElem, elementID,pProgMesh);
            break;

        case(ElementType::Beam):
            vProgElem.reserve(2);//分割された新しい要素の生成
            for(i=0; i< 2; i++){
                pProgElem= new CBeam; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();

                vProgElem.push_back(pProgElem);
            };
            dividBeam(pElem,vProgElem, elementID,pProgMesh);
            break;

        case(ElementType::Beam2):
            vProgElem.reserve(2);//分割された新しい要素の生成
            for(i=0; i< 2; i++){
                pProgElem= new CBeam2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                
                vProgElem.push_back(pProgElem);
            };
            dividBeam(pElem,vProgElem, elementID,pProgMesh);
            break;
    }//switch エンド
}

// 要素を分割して生成, for refineMesh()
// --
// Hexa(6面体)の分割
//
void CMeshFactory::dividHexa(CElement* pElem, vector<CElement*>& vProgElem, uint& elementID, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;//頂点のノード
    vector<CNode*> vEdgeNode;//辺のノード
    vector<CNode*> vFaceNode;//面のノード
    CNode          *pVolNode;//体中心のノード
    
    //    uint numOfVert,numOfEdge,numOfFace;//各要素の属性(分割の際に使用)
    //    numOfVert= NumberOfVertex::Hexa(); numOfFace= NumberOfFace::Hexa(); numOfEdge= NumberOfEdge::Hexa();

    uint i;
    //頂点のノード
    vVertNode.resize(8); for(i=0; i< 8; i++){ vVertNode[i] = pElem->getNode(i);}
    //辺のノード
    vEdgeNode.resize(12);for(i=0; i< 12; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    //面のノード
    vFaceNode.resize(6); for(i=0; i< 6; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    //体ノード
    pVolNode = pElem->getVolumeNode();


    // 8個のHexaを生成
    // 要素 0
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    vProgElem[0]->setNode(vFaceNode[0],2); vProgElem[0]->setNode(vEdgeNode[3],3);
    vProgElem[0]->setNode(vEdgeNode[8],4); vProgElem[0]->setNode(vFaceNode[4],5);
    vProgElem[0]->setNode(pVolNode,    6); vProgElem[0]->setNode(vFaceNode[3],7);

    pElem->setProgElem(vProgElem[0], 0);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット

    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    vProgElem[1]->setNode(vEdgeNode[1],2); vProgElem[1]->setNode(vFaceNode[0],3);
    vProgElem[1]->setNode(vFaceNode[4],4); vProgElem[1]->setNode(vEdgeNode[9],5);
    vProgElem[1]->setNode(vFaceNode[2],6); vProgElem[1]->setNode(pVolNode,    7);

    pElem->setProgElem(vProgElem[1], 1);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット

    // 要素 2
    vProgElem[2]->setNode(vEdgeNode[8],0); vProgElem[2]->setNode(vFaceNode[4],1);
    vProgElem[2]->setNode(pVolNode,    2); vProgElem[2]->setNode(vFaceNode[3],3);
    vProgElem[2]->setNode(vVertNode[4],4); vProgElem[2]->setNode(vEdgeNode[4],5);
    vProgElem[2]->setNode(vFaceNode[1],6); vProgElem[2]->setNode(vEdgeNode[7],7);

    pElem->setProgElem(vProgElem[2], 4);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット

    // 要素 3
    vProgElem[3]->setNode(vFaceNode[4],0); vProgElem[3]->setNode(vEdgeNode[9],1);
    vProgElem[3]->setNode(vFaceNode[2],2); vProgElem[3]->setNode(pVolNode,    3);
    vProgElem[3]->setNode(vEdgeNode[4],4); vProgElem[3]->setNode(vVertNode[5],5);
    vProgElem[3]->setNode(vEdgeNode[5],6); vProgElem[3]->setNode(vFaceNode[1],7);
    
    pElem->setProgElem(vProgElem[3], 5);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット

    // 要素 4
    vProgElem[4]->setNode(vEdgeNode[3],0); vProgElem[4]->setNode(vFaceNode[0],1);
    vProgElem[4]->setNode(vEdgeNode[2],2); vProgElem[4]->setNode(vVertNode[3],3);
    vProgElem[4]->setNode(vFaceNode[3],4); vProgElem[4]->setNode(pVolNode,    5);
    vProgElem[4]->setNode(vFaceNode[5],6); vProgElem[4]->setNode(vEdgeNode[11],7);
    
    pElem->setProgElem(vProgElem[4], 3);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット

    // 要素 5
    vProgElem[5]->setNode(vFaceNode[0],0); vProgElem[5]->setNode(vEdgeNode[1],1);
    vProgElem[5]->setNode(vVertNode[2],2); vProgElem[5]->setNode(vEdgeNode[2],3);
    vProgElem[5]->setNode(pVolNode,    4); vProgElem[5]->setNode(vFaceNode[2],5);
    vProgElem[5]->setNode(vEdgeNode[10],6);vProgElem[5]->setNode(vFaceNode[5],7);
    
    pElem->setProgElem(vProgElem[5], 2);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット

    // 要素 6
    vProgElem[6]->setNode(vFaceNode[3],0); vProgElem[6]->setNode(pVolNode,    1);
    vProgElem[6]->setNode(vFaceNode[5],2); vProgElem[6]->setNode(vEdgeNode[11],3);
    vProgElem[6]->setNode(vEdgeNode[7],4); vProgElem[6]->setNode(vFaceNode[1],5);
    vProgElem[6]->setNode(vEdgeNode[6],6); vProgElem[6]->setNode(vVertNode[7],7);
    
    pElem->setProgElem(vProgElem[6], 7);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット

    // 要素 7
    vProgElem[7]->setNode(pVolNode,    0); vProgElem[7]->setNode(vFaceNode[2],1);
    vProgElem[7]->setNode(vEdgeNode[10],2); vProgElem[7]->setNode(vFaceNode[5],3);
    vProgElem[7]->setNode(vFaceNode[1],4); vProgElem[7]->setNode(vEdgeNode[5],5);
    vProgElem[7]->setNode(vVertNode[6],6); vProgElem[7]->setNode(vEdgeNode[6],7);
    
    pElem->setProgElem(vProgElem[7], 6);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット

    
    // IDのセット
    for(i=0; i< 8; i++){
        //vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット
        vProgElem[i]->setID(elementID);          // <= elementIDは,いままでの数を数えているので,直前の配列数
        ++elementID;//直前の番号を渡した後でelementIDをカウントアップ

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iface;
    //マスター面
    if(pElem->isMPCMaster()){
        for(iface=0; iface< 6; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        setProgHexaMPCMaster(vProgElem, iface, 0,1,4,5);//Face_0 に生成される子要素:vProgElemの配列番号0,1,4,5
                        break;
                    case(1):
                        setProgHexaMPCMaster(vProgElem, iface, 2,3,6,7);//Face_1 に生成される子要素:vProgElemの配列番号2,3,6,7
                        break;
                    case(2):
                        setProgHexaMPCMaster(vProgElem, iface, 1,3,5,7);//Face_2 に生成される子要素:vProgElemの配列番号1,3,5,7
                        break;
                    case(3):
                        setProgHexaMPCMaster(vProgElem, iface, 0,2,4,6);//Face_3 に生成される子要素:vProgElemの配列番号0,2,4,6
                        break;
                    case(4):
                        setProgHexaMPCMaster(vProgElem, iface, 0,1,2,3);//Face_4 に生成される子要素:vProgElemの配列番号0,1,2,3
                        break;
                    case(5):
                        setProgHexaMPCMaster(vProgElem, iface, 4,5,6,7);//Face_5 に生成される子要素:vProgElemの配列番号4,5,6,7
                        break;
                }
            }
        };
    }
    //スレーブ面
    if(pElem->isMPCSlave()){
        for(iface=0; iface< 6; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        setProgHexaMPCSlave(vProgElem, iface, 0,1,4,5);
                        break;
                    case(1):
                        setProgHexaMPCSlave(vProgElem, iface, 2,3,6,7);
                        break;
                    case(2):
                        setProgHexaMPCSlave(vProgElem, iface, 1,3,5,7);
                        break;
                    case(3):
                        setProgHexaMPCSlave(vProgElem, iface, 0,2,4,6);
                        break;
                    case(4):
                        setProgHexaMPCSlave(vProgElem, iface, 0,1,2,3);
                        break;
                    case(5):
                        setProgHexaMPCSlave(vProgElem, iface, 4,5,6,7);
                        break;
                }
            }
        };
    }
    
    //CommMesh2(節点通信界面)
    if(pElem->isCommMesh2()){
        for(iface=0; iface< 6; iface++){
            if(pElem->isCommEntity(iface)){
                switch(iface){
                    case(0):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 0,1,4,5);//Face_0 に生成される子要素:vProgElemの配列番号0,1,4,5
                        break;
                    case(1):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 2,3,6,7);//Face_1 に生成される子要素:vProgElemの配列番号2,3,6,7
                        break;
                    case(2):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 1,3,5,7);//Face_2 に生成される子要素:vProgElemの配列番号1,3,5,7
                        break;
                    case(3):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 0,2,4,6);//Face_3 に生成される子要素:vProgElemの配列番号0,2,4,6
                        break;
                    case(4):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 0,1,2,3);//Face_4 に生成される子要素:vProgElemの配列番号0,1,2,3
                        break;
                    case(5):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 4,5,6,7);//Face_5 に生成される子要素:vProgElemの配列番号4,5,6,7
                        break;
                }
            }
        };
    }
    
}

// MPC面の属性をprogElemにセット
// --
// マスター
void CMeshFactory::setProgHexaMPCMaster(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l)
{
    vProgElem[i]->markingMPCMaster(); vProgElem[i]->markingMPCFace(iface);
    vProgElem[j]->markingMPCMaster(); vProgElem[j]->markingMPCFace(iface);
    vProgElem[k]->markingMPCMaster(); vProgElem[k]->markingMPCFace(iface);
    vProgElem[l]->markingMPCMaster(); vProgElem[l]->markingMPCFace(iface);
}
// スレーブ
void CMeshFactory::setProgHexaMPCSlave(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l)
{
    vProgElem[i]->markingMPCSlave(); vProgElem[i]->markingMPCFace(iface);
    vProgElem[j]->markingMPCSlave(); vProgElem[j]->markingMPCFace(iface);
    vProgElem[k]->markingMPCSlave(); vProgElem[k]->markingMPCFace(iface);
    vProgElem[l]->markingMPCSlave(); vProgElem[l]->markingMPCFace(iface);
}

// CommMesh2(節点通信界面)の属性をprogElemにセット
// --
void CMeshFactory::setProgHexaCommMesh2Entity(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l)
{
    vProgElem[i]->markingCommMesh2(); vProgElem[i]->markingCommEntity(iface);
    vProgElem[j]->markingCommMesh2(); vProgElem[j]->markingCommEntity(iface);
    vProgElem[k]->markingCommMesh2(); vProgElem[k]->markingCommEntity(iface);
    vProgElem[l]->markingCommMesh2(); vProgElem[l]->markingCommEntity(iface);
}

// 4面体の分割
//
void CMeshFactory::dividTetra(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;//頂点のノード
    vector<CNode*> vEdgeNode;//辺のノード
    vector<CNode*> vFaceNode;//面のノード
    CNode          *pVolNode;//体中心のノード
    
    uint numOfVert,numOfEdge,numOfFace;//各要素の属性(分割の際に使用)
    numOfVert= NumberOfVertex::Tetra(); numOfFace= NumberOfFace::Tetra(); numOfEdge= NumberOfEdge::Tetra();

    uint i;
    //頂点のノード
    vVertNode.resize(numOfVert);
    for(i=0; i< numOfVert; i++){
        vVertNode[i] = pElem->getNode(i);
    }
    //辺のノード
    vEdgeNode.resize(numOfEdge);
    for(i=0; i< numOfEdge; i++){
        vEdgeNode[i] = pElem->getEdgeInterNode(i);
    }
    //面のノード
    vFaceNode.resize(numOfFace);
    for(i=0; i< numOfFace; i++){
        vFaceNode[i] = pElem->getFaceNode(i);
    }
    //体ノード
    pVolNode = pElem->getVolumeNode();


    // 4個のHexaを生成
    // --
    // 要素 0
    vProgElem[0]->setNode(vEdgeNode[2],0); vProgElem[0]->setNode(vVertNode[0],1);
    vProgElem[0]->setNode(vEdgeNode[0],2); vProgElem[0]->setNode(vFaceNode[0],3);
    vProgElem[0]->setNode(vFaceNode[3],4); vProgElem[0]->setNode(vEdgeNode[3],5);
    vProgElem[0]->setNode(vFaceNode[1],6); vProgElem[0]->setNode(pVolNode,    7);
    
    pElem->setProgElem(vProgElem[0], 0);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット
    
    // 要素 1
    vProgElem[1]->setNode(vVertNode[2],0); vProgElem[1]->setNode(vEdgeNode[2],1);
    vProgElem[1]->setNode(vFaceNode[0],2); vProgElem[1]->setNode(vEdgeNode[1],3);
    vProgElem[1]->setNode(vEdgeNode[5],4); vProgElem[1]->setNode(vFaceNode[3],5);
    vProgElem[1]->setNode(pVolNode,    6); vProgElem[1]->setNode(vFaceNode[2],7);
    
    pElem->setProgElem(vProgElem[1], 2);
    
    // 要素 2
    vProgElem[2]->setNode(vFaceNode[0],0); vProgElem[2]->setNode(vEdgeNode[0],1);
    vProgElem[2]->setNode(vVertNode[1],2); vProgElem[2]->setNode(vEdgeNode[1],3);
    vProgElem[2]->setNode(pVolNode,    4); vProgElem[2]->setNode(vFaceNode[1],5);
    vProgElem[2]->setNode(vEdgeNode[4],6); vProgElem[2]->setNode(vFaceNode[2],7);
    
    pElem->setProgElem(vProgElem[2], 1);
    
    // 要素 3
    vProgElem[3]->setNode(vFaceNode[3],0); vProgElem[3]->setNode(vEdgeNode[3],1);
    vProgElem[3]->setNode(vFaceNode[1],2); vProgElem[3]->setNode(pVolNode,    3);
    vProgElem[3]->setNode(vEdgeNode[5],4); vProgElem[3]->setNode(vVertNode[3],5);
    vProgElem[3]->setNode(vEdgeNode[4],6); vProgElem[3]->setNode(vFaceNode[2],7);
    
    pElem->setProgElem(vProgElem[3], 3);

    ////debug
    //cout << "要素にノードをセット@dividTetra" << endl;

    // IDのセット
    for(i=0; i< 4; i++){
        //vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iface;
    //マスター面
    if(pElem->isMPCMaster()){
        for(iface=0; iface< 4; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(0);
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(0);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(0);
                        break;
                    case(1):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(2);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(2);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(2);
                        break;
                    case(2):
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(3);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(5);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(1);
                        break;
                    case(3):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(4);
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(4);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(4);
                        break;
                }
            }
        };
    }
    //スレーブ面
    if(pElem->isMPCSlave()){
        for(iface=0; iface< 4; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(0);
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(0);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(0);
                        break;
                    case(1):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(2);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(2);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(2);
                        break;
                    case(2):
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(3);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(5);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(1);
                        break;
                    case(3):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(4);
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(4);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(4);
                        break;
                }
            }
        };
    }

    //CommMesh2(節点通信界面)
    if(pElem->isCommMesh2()){
        for(iface=0; iface< 4; iface++){
            if(pElem->isCommEntity(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(0);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(0);
                        break;
                    case(1):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(2);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(2);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(2);
                        break;
                    case(2):
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(3);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(5);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(1);
                        break;
                    case(3):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(4);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(4);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(4);
                        break;
                }
            }
        };
    }
}


// プリズムの分割
//
void CMeshFactory::dividPrism(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;//頂点のノード
    vector<CNode*> vEdgeNode;//辺のノード
    vector<CNode*> vFaceNode;//面のノード
    CNode          *pVolNode;//体中心のノード
    
    uint numOfVert,numOfEdge,numOfFace;//各要素の属性(分割の際に使用)
    numOfVert= NumberOfVertex::Prism(); numOfFace= NumberOfFace::Prism(); numOfEdge= NumberOfEdge::Prism();

    uint i;
    //頂点のノード
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    //辺のノード
    vEdgeNode.resize(numOfEdge);for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    //面のノード
    vFaceNode.resize(numOfFace); for(i=0; i< numOfFace; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    //体ノード
    pVolNode = pElem->getVolumeNode();
    
    // 要素 0
    vProgElem[0]->setNode(vVertNode[2],0); vProgElem[0]->setNode(vEdgeNode[1],1);
    vProgElem[0]->setNode(vFaceNode[0],2); vProgElem[0]->setNode(vEdgeNode[2],3);
    vProgElem[0]->setNode(vEdgeNode[5],4); vProgElem[0]->setNode(vFaceNode[4],5);
    vProgElem[0]->setNode(pVolNode,    6); vProgElem[0]->setNode(vFaceNode[3],7);
    
    pElem->setProgElem(vProgElem[0], 2);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット
    
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[1],0); vProgElem[1]->setNode(vVertNode[0],1);
    vProgElem[1]->setNode(vEdgeNode[0],2); vProgElem[1]->setNode(vFaceNode[0],3);
    vProgElem[1]->setNode(vFaceNode[4],4); vProgElem[1]->setNode(vEdgeNode[3],5);
    vProgElem[1]->setNode(vFaceNode[2],6); vProgElem[1]->setNode(pVolNode,    7);
    
    pElem->setProgElem(vProgElem[1], 0);
    
    // 要素 2
    vProgElem[2]->setNode(vFaceNode[0],0); vProgElem[2]->setNode(vEdgeNode[0],1);
    vProgElem[2]->setNode(vVertNode[1],2); vProgElem[2]->setNode(vEdgeNode[2],3);
    vProgElem[2]->setNode(pVolNode,    4); vProgElem[2]->setNode(vFaceNode[2],5);
    vProgElem[2]->setNode(vEdgeNode[4],6); vProgElem[2]->setNode(vFaceNode[3],7);
    
    pElem->setProgElem(vProgElem[2], 1);

    // 要素 3
    vProgElem[3]->setNode(vEdgeNode[5],0); vProgElem[3]->setNode(vFaceNode[4],1);
    vProgElem[3]->setNode(pVolNode,    2); vProgElem[3]->setNode(vFaceNode[3],3);
    vProgElem[3]->setNode(vVertNode[5],4); vProgElem[3]->setNode(vEdgeNode[8],5);
    vProgElem[3]->setNode(vFaceNode[1],6); vProgElem[3]->setNode(vEdgeNode[7],7);
    
    pElem->setProgElem(vProgElem[3], 5);
    
    // 要素 4
    vProgElem[4]->setNode(vFaceNode[4],0); vProgElem[4]->setNode(vEdgeNode[3],1);
    vProgElem[4]->setNode(vFaceNode[2],2); vProgElem[4]->setNode(pVolNode,    3);
    vProgElem[4]->setNode(vEdgeNode[8],4); vProgElem[4]->setNode(vVertNode[3],5);
    vProgElem[4]->setNode(vEdgeNode[6],6); vProgElem[4]->setNode(vFaceNode[1],7);
    
    pElem->setProgElem(vProgElem[4], 3);
    
    // 要素 5
    vProgElem[5]->setNode(pVolNode,    0); vProgElem[5]->setNode(vFaceNode[2],1);
    vProgElem[5]->setNode(vEdgeNode[4],2); vProgElem[5]->setNode(vFaceNode[3],3);
    vProgElem[5]->setNode(vFaceNode[1],4); vProgElem[5]->setNode(vEdgeNode[6],5);
    vProgElem[5]->setNode(vVertNode[4],6); vProgElem[5]->setNode(vEdgeNode[7],7);
    
    pElem->setProgElem(vProgElem[5], 4);

    // IDのセット
    for(i=0; i< 6; i++){
        //vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iface;
    //マスター面
    if(pElem->isMPCMaster()){
        for(iface=0; iface< 5; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(0);
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(0);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(0);
                        break;
                    case(1):
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(1);
                        vProgElem[4]->markingMPCMaster(); vProgElem[4]->markingMPCFace(1);
                        vProgElem[5]->markingMPCMaster(); vProgElem[5]->markingMPCFace(1);
                        break;
                    case(2):
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(2);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(2);
                        vProgElem[4]->markingMPCMaster(); vProgElem[4]->markingMPCFace(2);
                        vProgElem[5]->markingMPCMaster(); vProgElem[5]->markingMPCFace(2);
                        break;
                    case(3):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(3);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(5);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(3);
                        vProgElem[5]->markingMPCMaster(); vProgElem[5]->markingMPCFace(5);
                        break;
                    case(4):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(4);
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(4);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(4);
                        vProgElem[4]->markingMPCMaster(); vProgElem[4]->markingMPCFace(4);
                        break;
                }
            }
        };
    }
    //スレーブ面
    if(pElem->isMPCSlave()){
        for(iface=0; iface< 5; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(0);
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(0);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(0);
                        break;
                    case(1):
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(1);
                        vProgElem[4]->markingMPCSlave(); vProgElem[4]->markingMPCFace(1);
                        vProgElem[5]->markingMPCSlave(); vProgElem[5]->markingMPCFace(1);
                        break;
                    case(2):
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(2);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(2);
                        vProgElem[4]->markingMPCSlave(); vProgElem[4]->markingMPCFace(2);
                        vProgElem[5]->markingMPCSlave(); vProgElem[5]->markingMPCFace(2);
                        break;
                    case(3):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(3);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(5);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(3);
                        vProgElem[5]->markingMPCSlave(); vProgElem[5]->markingMPCFace(5);
                        break;
                    case(4):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(4);
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(4);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(4);
                        vProgElem[4]->markingMPCSlave(); vProgElem[4]->markingMPCFace(4);
                        break;
                }
            }
        };
    }

    //CommMesh2(節点通信界面)
    if(pElem->isCommMesh2()){
        for(iface=0; iface< 5; iface++){
            if(pElem->isCommEntity(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(0);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(0);
                        break;
                    case(1):
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(1);
                        vProgElem[4]->markingCommMesh2(); vProgElem[4]->markingCommEntity(1);
                        vProgElem[5]->markingCommMesh2(); vProgElem[5]->markingCommEntity(1);
                        break;
                    case(2):
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(2);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(2);
                        vProgElem[4]->markingCommMesh2(); vProgElem[4]->markingCommEntity(2);
                        vProgElem[5]->markingCommMesh2(); vProgElem[5]->markingCommEntity(2);
                        break;
                    case(3):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(3);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(5);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(3);
                        vProgElem[5]->markingCommMesh2(); vProgElem[5]->markingCommEntity(5);
                        break;
                    case(4):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(4);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(4);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(4);
                        vProgElem[4]->markingCommMesh2(); vProgElem[4]->markingCommEntity(4);
                        break;
                }
            }
        };
    }
}
// ピラミッドの分割 <= 削除
//
void CMeshFactory::dividPyramid(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;//頂点のノード
    vector<CNode*> vEdgeNode;//辺のノード
    vector<CNode*> vFaceNode;//面のノード
    CNode          *pVolNode;//体中心のノード
    
    uint numOfVert,numOfEdge,numOfFace;//各要素の属性(分割の際に使用)
    numOfVert= NumberOfVertex::Pyramid(); numOfFace= NumberOfFace::Pyramid(); numOfEdge= NumberOfEdge::Pyramid();

    uint i;
    //頂点のノード
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    //辺のノード
    vEdgeNode.resize(numOfEdge);for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    //面のノード
    vFaceNode.resize(numOfFace); for(i=0; i< numOfFace; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    //体ノード
    pVolNode = pElem->getVolumeNode();
    
    // 要素 0
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    vProgElem[0]->setNode(vFaceNode[0],2); vProgElem[0]->setNode(vEdgeNode[3],3);
    vProgElem[0]->setNode(vEdgeNode[7],4); vProgElem[0]->setNode(vFaceNode[4],5);
    vProgElem[0]->setNode(pVolNode,    6); vProgElem[0]->setNode(vFaceNode[3],7);
    
    pElem->setProgElem(vProgElem[0], 0);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット
    
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    vProgElem[1]->setNode(vEdgeNode[1],2); vProgElem[1]->setNode(vFaceNode[0],3);
    vProgElem[1]->setNode(vFaceNode[4],4); vProgElem[1]->setNode(vEdgeNode[4],5);
    vProgElem[1]->setNode(vFaceNode[1],6); vProgElem[1]->setNode(pVolNode,    7);
    
    pElem->setProgElem(vProgElem[1], 1);
    
    // 要素 2
    vProgElem[2]->setNode(vFaceNode[0],0); vProgElem[2]->setNode(vEdgeNode[1],1);
    vProgElem[2]->setNode(vVertNode[2],2); vProgElem[2]->setNode(vEdgeNode[2],3);
    vProgElem[2]->setNode(pVolNode,    4); vProgElem[2]->setNode(vFaceNode[1],5);
    vProgElem[2]->setNode(vEdgeNode[5],6); vProgElem[2]->setNode(vFaceNode[2],7);
    
    pElem->setProgElem(vProgElem[2], 2);
    
    // 要素 3
    vProgElem[3]->setNode(vEdgeNode[3],0); vProgElem[3]->setNode(vFaceNode[0],1);
    vProgElem[3]->setNode(vEdgeNode[2],2); vProgElem[3]->setNode(vVertNode[3],3);
    vProgElem[3]->setNode(vFaceNode[3],4); vProgElem[3]->setNode(pVolNode,    5);
    vProgElem[3]->setNode(vFaceNode[2],6); vProgElem[3]->setNode(vEdgeNode[6],7);
    
    pElem->setProgElem(vProgElem[3], 3);
    

//    // 要素 4 (pyramid)
//    vProgElem[4]->setNode(vEdgeNode[7],0); vProgElem[4]->setNode(vFaceNode[4],1);
//    vProgElem[4]->setNode(pVolNode,    2); vProgElem[4]->setNode(vFaceNode[3],3);
//    vProgElem[4]->setNode(vVertNode[4],4);
//    // 要素 5 (pyramid)
//    vProgElem[5]->setNode(vFaceNode[4],0); vProgElem[5]->setNode(vEdgeNode[4],1);
//    vProgElem[5]->setNode(vFaceNode[1],2); vProgElem[5]->setNode(pVolNode,    3);
//    vProgElem[5]->setNode(vVertNode[4],4);
//    // 要素 6 (pyramid)
//    vProgElem[6]->setNode(vFaceNode[1],0); vProgElem[6]->setNode(vEdgeNode[5],1);
//    vProgElem[6]->setNode(vFaceNode[2],2); vProgElem[6]->setNode(pVolNode,    3);
//    vProgElem[6]->setNode(vVertNode[4],4);
//    // 要素 7 (pyramid)
//    vProgElem[7]->setNode(vFaceNode[2],0); vProgElem[7]->setNode(vEdgeNode[6],1);
//    vProgElem[7]->setNode(vFaceNode[3],2); vProgElem[7]->setNode(pVolNode,    3);
//    vProgElem[7]->setNode(vVertNode[4],4);


    // 要素 4 (pyramid)
    vProgElem[4]->setNode(vEdgeNode[7],0); vProgElem[4]->setNode(vVertNode[4],1);
    vProgElem[4]->setNode(vEdgeNode[4],2); vProgElem[4]->setNode(vFaceNode[4],3);
    vProgElem[4]->setNode(pVolNode,    4);
    
    // pyramid分割されたProgElemは全て頂点4を所有するので,
    // 頂点番号によるProgElemの分類に加えて,Face順に数える.
    pElem->setProgElem(vProgElem[4], 7);//FaceNode[4]
    
    // 要素 5 (pyramid)
    vProgElem[5]->setNode(vEdgeNode[4],0); vProgElem[5]->setNode(vVertNode[4],1);
    vProgElem[5]->setNode(vEdgeNode[5],2); vProgElem[5]->setNode(vFaceNode[1],3);
    vProgElem[5]->setNode(pVolNode,    4);
    
    pElem->setProgElem(vProgElem[5], 4);//FaceNode[1]
    
    // 要素 6 (pyramid)
    vProgElem[6]->setNode(vEdgeNode[5],0); vProgElem[6]->setNode(vVertNode[4],1);
    vProgElem[6]->setNode(vEdgeNode[6],2); vProgElem[6]->setNode(vFaceNode[2],3);
    vProgElem[6]->setNode(pVolNode,    4);
    
    pElem->setProgElem(vProgElem[6], 5);//FaceNode[2]
    
    // 要素 7 (pyramid)
    vProgElem[7]->setNode(vEdgeNode[6],0); vProgElem[7]->setNode(vVertNode[4],1);
    vProgElem[7]->setNode(vEdgeNode[7],2); vProgElem[7]->setNode(vFaceNode[3],3);
    vProgElem[7]->setNode(pVolNode,    4);
    
    pElem->setProgElem(vProgElem[7], 6);//FaceNode[3]

    // IDのセット
    for(i=0; i< 8; i++){
        //vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット
        
        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    //MPCは実装しない=> ピラミッドは削除予定

}
// 四辺形の分割
//
void CMeshFactory::dividQuad(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;//頂点のノード
    vector<CNode*> vEdgeNode;//辺のノード
    vector<CNode*> vFaceNode;//面のノード
    
    uint numOfVert,numOfEdge,numOfFace;//各要素の属性(分割の際に使用)
    numOfVert= NumberOfVertex::Quad(); numOfFace= NumberOfFace::Quad(); numOfEdge= NumberOfEdge::Quad();

    uint i;
    //頂点のノード
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    //辺のノード
    vEdgeNode.resize(numOfEdge);for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    //面のノード
    vFaceNode.resize(numOfFace); for(i=0; i< numOfFace; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    
    // 要素 0
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    vProgElem[0]->setNode(vFaceNode[0],2); vProgElem[0]->setNode(vEdgeNode[3],3);
    
    pElem->setProgElem(vProgElem[0], 0);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット
    
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    vProgElem[1]->setNode(vEdgeNode[1],2); vProgElem[1]->setNode(vFaceNode[0],3);
    
    pElem->setProgElem(vProgElem[1], 1);
    
    // 要素 2
    vProgElem[2]->setNode(vEdgeNode[1],0); vProgElem[2]->setNode(vVertNode[2],1);
    vProgElem[2]->setNode(vEdgeNode[2],2); vProgElem[2]->setNode(vFaceNode[0],3);
    
    pElem->setProgElem(vProgElem[2], 2);
    
    // 要素 3
    vProgElem[3]->setNode(vEdgeNode[2],0); vProgElem[3]->setNode(vVertNode[3],1);
    vProgElem[3]->setNode(vEdgeNode[3],2); vProgElem[3]->setNode(vFaceNode[0],3);
    
    pElem->setProgElem(vProgElem[3], 3);

    // IDのセット
    for(i=0; i< 4; i++){
        //vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iprog;
    //マスター面
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 4; iprog++){ vProgElem[iprog]->markingMPCMaster(); vProgElem[iprog]->markingMPCFace(0);}
    }
    //スレーブ面
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 4; iprog++){ vProgElem[iprog]->markingMPCSlave(); vProgElem[iprog]->markingMPCFace(0);}
    }
    
    //CommMesh2(節点通信界面):Quadなので通信"辺"
    uint iedge;
    if(pElem->isCommMesh2()){
        for(iedge=0; iedge< 4; iedge++){
            if(pElem->isCommEntity(iedge)){
                switch(iedge){
                    case(0):
                        //辺番号"0" にくっついているprogElemは,progElem[0]とprogElem[1]
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);//prog[0]の辺=0
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(0);//prog[1]の辺=0
                        break;
                    case(1):
                        //辺番号"1" にくっついているproElemは,progElem[1]とprogElem[2]
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(1);//prog[1]の辺=1
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(0);//prog[2]の辺=0
                        break;
                    case(2):
                        //辺番号"2" にくっついているproElemは,progElem[2]とprogElem[3]
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(1);//prog[2]の辺=1
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(0);//prog[3]の辺=0
                        break;
                    case(3):
                        //辺番号"3" にくっついているproElemは,progElem[3]とprogElem[0]
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(1);//prog[3]の辺=1
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(3);//prog[0]の辺=3
                        break;
                }
            }//if( CommEntity )
        };//for(iedge)
    }//if(CommMesh2)
    
}
// 三角形の分割
//
void CMeshFactory::dividTriangle(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;//頂点のノード
    vector<CNode*> vEdgeNode;//辺のノード
    vector<CNode*> vFaceNode;//面のノード
    
    uint numOfVert,numOfEdge,numOfFace;//各要素の属性(分割の際に使用)
    numOfVert= NumberOfVertex::Triangle(); numOfFace= NumberOfFace::Triangle(); numOfEdge= NumberOfEdge::Triangle();

    uint i;
    //頂点のノード
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    //辺のノード
    vEdgeNode.resize(numOfEdge);for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    //面のノード
    vFaceNode.resize(numOfFace); for(i=0; i< numOfFace; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    
    // 要素 0
    vProgElem[0]->setNode(vEdgeNode[0],0); vProgElem[0]->setNode(vVertNode[1],1);
    vProgElem[0]->setNode(vEdgeNode[1],2); vProgElem[0]->setNode(vFaceNode[0],3);
    
    pElem->setProgElem(vProgElem[0], 1);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット
    
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[1],0); vProgElem[1]->setNode(vVertNode[2],1);
    vProgElem[1]->setNode(vEdgeNode[2],2); vProgElem[1]->setNode(vFaceNode[0],3);
    
    pElem->setProgElem(vProgElem[1], 2);
    
    // 要素 2
    vProgElem[2]->setNode(vEdgeNode[2],0); vProgElem[2]->setNode(vVertNode[0],1);
    vProgElem[2]->setNode(vEdgeNode[0],2); vProgElem[2]->setNode(vFaceNode[0],3);
    
    pElem->setProgElem(vProgElem[2], 0);

    // IDのセット
    for(i=0; i< 3; i++){
        //vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iprog;
    // マスター面
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 3; iprog++){ vProgElem[iprog]->markingMPCMaster(); vProgElem[iprog]->markingMPCFace(0);}
    }
    // スレーブ面
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 3; iprog++){ vProgElem[iprog]->markingMPCSlave(); vProgElem[iprog]->markingMPCFace(0);}
    }

    // CommMesh2(通信節点界面) : Triなので通信"辺"
    uint iedge;
    if(pElem->isCommMesh2()){
        for(iedge=0; iedge< 3; iedge++){
            if(pElem->isCommEntity(iedge)){
                switch(iedge){
                    case(0):
                        //辺"0"に生成されるprogElemは,prog[2]とprog[0]
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(1);//progElemの辺”1”
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);//progElemの辺”0”
                        break;
                    case(1):
                        //辺"1"に生成されるprogElemは,prog[0]とprog[1]
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(1);//progElemの辺”1”
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(0);//progElemの辺”0”
                        break;
                    case(2):
                        //辺"2"に生成されるprogElemは,prog[1]とprog[2]
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(1);//progElemの辺”1”
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(0);//progElemの辺”0”
                        break;
                }
            }// if( CommEntity )
        };
    }// if(CommMesh2)
}
// ビームの分割
//
void CMeshFactory::dividBeam(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;//頂点のノード
    vector<CNode*> vEdgeNode;//辺のノード
    
    uint numOfVert,numOfEdge;//各要素の属性(分割の際に使用)
    numOfVert= NumberOfVertex::Beam(); numOfEdge= NumberOfEdge::Beam();

    uint i;
    //頂点のノード
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    //辺のノード
    vEdgeNode.resize(numOfEdge); for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    
    // 要素 0
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    
    pElem->setProgElem(vProgElem[0], 0);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット
    
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);

    pElem->setProgElem(vProgElem[1], 1);//CommElemのために頂点番号(vVertNode)順に親ElemにProgElemをセット

    // IDのセット
    for(i=0; i< 2; i++){
        //vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iprog;
    // マスター
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 2; iprog++) vProgElem[iprog]->markingMPCMaster();
    }
    // スレーブ
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 2; iprog++) vProgElem[iprog]->markingMPCSlave();
    }

    // CommMesh2(節点通信界面) :Beamなので通信”点” { progしても変化なし.}
    uint ivert;
    if(pElem->isCommMesh2()){
        for(ivert=0; ivert< 2; ivert++){
            if(pElem->isCommEntity(ivert)){
                switch(ivert){
                    case(0):
                        //頂点"0"に生成されるprogElem, progElemの頂点"0"が通信点
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);
                        break;
                    case(1):
                        //頂点"1"に生成されるprogElem, progElemの頂点"1"が通信点
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(1);
                        break;
                }
            }//if(CommEntity)
        };
    }// if(CommMesh2)
}



// setup to BucketMesh in AssyModel
//
void CMeshFactory::setupBucketMesh(const uint& mgLevel, const uint& maxID, const uint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);

    mpTAssyModel->intializeBucket(maxID, minID);

    mpTAssyModel->setMaxMeshID(maxID);
    mpTAssyModel->setMinMeshID(minID);
}


// GMGModel MultiGridの各層のAssyModelを全て生成(AssyModel変数を確保)
//
void CMeshFactory::GeneAssyModel(const uint& num_of_mgLevel)
{
    mpGMGModel->initAssyModel(num_of_mgLevel);

    mpGMGModel->reserveAssyModel(num_of_mgLevel);
    
    for(uint i=0; i<num_of_mgLevel; i++){
        mpTAssyModel = new CAssyModel();
        mpTAssyModel->setMGLevel(i);

        mpGMGModel->addModel(mpTAssyModel,i);
    };
}

// Mesh リザーブ
//
void CMeshFactory::reserveMesh(const uint& mgLevel, const uint& num_of_mesh)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTAssyModel->resizeMesh(num_of_mesh);
}

// Mesh set to AssyModel
//
void CMeshFactory::GeneMesh(const uint& mgLevel, const uint& mesh_id, const uint& index, const uint& nProp)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);

    //AssyModelにMeshを生成してセット
    //----
    mpTMesh = new CMesh();

    mpTMesh->setMeshID(mesh_id);
    mpTMesh->setMGLevel(mgLevel);
    mpTMesh->setMaxMGLevel(mMGLevel);//Factoryのレベル==最大階層数
    mpTMesh->setSolutionType(mnSolutionType);
    mpTMesh->setProp(nProp);

    mpTAssyModel->setBucket(mesh_id, index);
    mpTAssyModel->setMesh(mpTMesh, index);

    //各MeshにBNodeMeshGrpを生成してセット
    //----
    CBNodeMeshGrp *pBNodeMeshGrp= new CBNodeMeshGrp;
    mpTMesh->setBNodeMeshGrp(pBNodeMeshGrp);
}


// Node Genetate 
//
void CMeshFactory::GeneNode(const uint& mgLevel, const uint& mesh_id, const uint& id, const vdouble& coord,
                            const uint& nodeType, const uint& numOfScaParam, const uint& numOfVecParam)
{
    CNode *pNode;

    switch(nodeType){
        case(NodeType::Scalar):
            pNode = new CScalarNode();
            pNode->resizeScalar(numOfScaParam);
            break;
        case(NodeType::Vector):
            pNode = new CVectorNode();
            pNode->resizeVector(numOfVecParam);
            break;
        case(NodeType::ScalarVector):
            pNode = new CScalarVectorNode();
            pNode->resizeScalar(numOfScaParam);
            pNode->resizeVector(numOfVecParam);
            break;
        default:
            //pNode->InitializeNodeADOF(vParam, num_of_param);
            break;
    }
    pNode->setID(id);// id は、連続したIndex番号にする予定(09.06.23)
    pNode->setCoord(coord);

    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);// MultiGrid Level==0

    if(!mpTAssyModel) mpLogger->Info(Utility::LoggerMode::MWDebug, "AssyModel => NULL");//debug

    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    ////debug
    //if(!mpTMesh) cout << "Factory::mpTMesh == NULL" << endl;

    //mpTMesh->setNode(pNode,id);
    mpTMesh->setNode(pNode);
}


// Mesh::seupNumOfNode
//
void CMeshFactory::setupNode(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    
    mpTMesh->setupNumOfNode();
}

// Mesh::reserveNode
//
void CMeshFactory::reserveNode(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->reserveNode(num_of_node);
}



// Element Generate
//
void CMeshFactory::GeneElement(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& type_num, const vint& node_id)
{
    CElement *pElement;

    switch(type_num){
        // 1次
        case(ElementType::Hexa):
            pElement = new CHexa;
            break;
        case(ElementType::Tetra):
            pElement = new CTetra;
            break;
        case(ElementType::Prism):
            pElement = new CPrism;
            break;
        case(ElementType::Quad):
            pElement = new CQuad;
            break;
        case(ElementType::Triangle):
            pElement = new CTriangle;
            break;
        case(ElementType::Beam):
            pElement = new CBeam;
            break;
        // 2次
        case(ElementType::Hexa2):
            pElement = new CHexa2;
            break;
        case(ElementType::Tetra2):
            pElement = new CTetra2;
            break;
        case(ElementType::Prism2):
            pElement = new CPrism2;
            break;
        case(ElementType::Quad2):
            pElement = new CQuad2;
            break;
        case(ElementType::Triangle2):
            pElement = new CTriangle2;
            break;
        case(ElementType::Beam2):
            pElement = new CBeam2;
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error,"Error::GeneElement at Factory");
            break;
    }
    pElement->initialize();
    pElement->setID(id);

    CNode *pNode;
    uint i;
    // Node* setup
    //
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);    
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    for(i=0; i < node_id.size(); i++){
        pNode = mpTMesh->getNode(node_id[i]);
        pElement->setNode(pNode, i);
    };
    //2次要素の場合"EdgeInterNode"にセット
    if(pElement->getOrder()==ElementOrder::Second){
        uint iedge;
        uint nNumOfVert = pElement->getNumOfVert();
        for(iedge=0; iedge < pElement->getNumOfEdge(); iedge++){

            uint nNodeID = node_id[nNumOfVert + iedge];

            pNode = mpTMesh->getNode(nNodeID);
            pElement->setEdgeInterNode(pNode, iedge);
        };
    }
    
    // Elem to Mesh
    //mpTMesh->setElement(pElement,id);
    mpTMesh->setElement(pElement);
}


// Mesh::setupNumOfElement
//
void CMeshFactory::setupElement(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupNumOfElement();
}

// Mesh::reserveElement
//
void CMeshFactory::reserveElement(const uint& mgLevel, const uint& mesh_id, const uint& num_of_element)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveElement(num_of_element);
}

// 節点まわりの要素集合
// Mesh::reserveAggElement
//
void CMeshFactory::resizeAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->resizeAggregate(num_of_node);
}

void CMeshFactory::GeneAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    uint i;
    if(mnSolutionType==SolutionType::FEM){
        for(i=0; i< num_of_node; i++){
            CAggregateElement *pAggElem = new CAggregateElement;

            mpTMesh->setAggElement(pAggElem, i);//Node周辺の要素集合
        };
    }
    if(mnSolutionType==SolutionType::FVM){
        for(i=0; i< num_of_node; i++){
            CAggregateNode    *pAggNode = new CAggregateNode;

            mpTMesh->setAggNode(pAggNode, i);//Node周辺のNode集合
        };
    }
}


//----
// Boundary
//----
void CMeshFactory::reserveBoundaryNodeMesh(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->reserveBndNodeMesh(num_of_bnd);
}
// BoundaryNodeMesh 生成
void CMeshFactory::GeneBoundaryNodeMesh(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id, const uint& bnd_type, const string& bnd_name)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    
    CBoundaryNodeMesh *pBNodeMesh= new CBoundaryNodeMesh;

    pBNodeMesh->setID(bnd_id);
    pBNodeMesh->setBndType(bnd_type);
    pBNodeMesh->setName(bnd_name);

    mpTMesh->setBndNodeMesh(pBNodeMesh);
}
uint CMeshFactory::getNumOfBounaryNodeMesh(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    return mpTMesh->getNumOfBoundaryNodeMesh();
}

void CMeshFactory::reserveBoundaryFaceMesh(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    
    mpTMesh->reserveBndFaceMesh(num_of_bnd);
}
// BoundaryFaceMesh 生成
void CMeshFactory::GeneBoundaryFaceMesh(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id, const uint& bnd_type, const string& bnd_name,
                                        const uint& numOfDOF, const vuint& vDOF)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    

    CBoundaryFaceMesh *pBFaceMesh= new CBoundaryFaceMesh;

    pBFaceMesh->setMGLevel(mgLevel);
    pBFaceMesh->setMaxMGLevel(mMGLevel);
    pBFaceMesh->setID(bnd_id);
    pBFaceMesh->setBndType(bnd_type);
    pBFaceMesh->setName(bnd_name);

    pBFaceMesh->resizeDOF(numOfDOF);

    if(numOfDOF != vDOF.size()) mpLogger->Info(Utility::LoggerMode::Error, "CMeshFactory::GeneBoundaryFaceMesh, invalid argument");

    uint idof, dof;
    for(idof=0; idof < vDOF.size(); idof++){
        dof = vDOF[idof];
        pBFaceMesh->setDOF(idof, dof);
    };

    mpTMesh->setBndFaceMesh(pBFaceMesh);

}
uint CMeshFactory::getNumOfBounaryFaceMesh(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    return mpTMesh->getNumOfBoundaryFaceMesh();
}

void CMeshFactory::reserveBoundaryVolumeMesh(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->reserveBndVolumeMesh(num_of_bnd);
}
// BoundaryVolumeMesh 生成
void CMeshFactory::GeneBoundaryVolumeMesh(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id, const uint& bnd_type, const string& bnd_name,
                                          const uint& numOfDOF, const vuint& vDOF)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);


    CBoundaryVolumeMesh *pBVolMesh= new CBoundaryVolumeMesh;
    
    pBVolMesh->setMGLevel(mgLevel);
    pBVolMesh->setMaxMGLevel(mMGLevel);
    pBVolMesh->setID(bnd_id);
    pBVolMesh->setBndType(bnd_type);
    pBVolMesh->setName(bnd_name);

    pBVolMesh->resizeDOF(numOfDOF);

    if(numOfDOF != vDOF.size()) mpLogger->Info(Utility::LoggerMode::Error, "CMeshFactory::GeneBoundaryVolumeMesh, invalid argument");

    uint idof, dof;
    for(idof=0; idof < vDOF.size(); idof++){
        dof = vDOF[idof];
        pBVolMesh->setDOF(idof, dof);
    }

    mpTMesh->setBndVolumeMesh(pBVolMesh);
    
}
uint CMeshFactory::getNumOfBounaryVolumeMesh(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    
    return mpTMesh->getNumOfBoundaryVolumeMesh();
}

void CMeshFactory::reserveBoundaryEdgeMesh(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->reserveBndEdgeMesh(num_of_bnd);
}
// BoundaryEdgeMesh 生成
void CMeshFactory::GeneBoundaryEdgeMesh(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id, const uint& bnd_type, const string& bnd_name,
                   const uint& numOfDOF, const vuint& vDOF)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    CBoundaryEdgeMesh *pBEdgeMesh= new CBoundaryEdgeMesh;
    
    pBEdgeMesh->setMGLevel(mgLevel);
    pBEdgeMesh->setMaxMGLevel(mMGLevel);
    pBEdgeMesh->setID(bnd_id);
    pBEdgeMesh->setBndType(bnd_type);
    pBEdgeMesh->setName(bnd_name);

    pBEdgeMesh->resizeDOF(numOfDOF);
    if(numOfDOF != vDOF.size()) mpLogger->Info(Utility::LoggerMode::Error, "CMeshFactory::GeneBoundaryEdgeMesh, invalid argument");

    uint idof, dof;
    for(idof=0; idof < numOfDOF; idof++){
        dof = vDOF[idof];
        pBEdgeMesh->setDOF(idof, dof);
    }

    mpTMesh->setBndEdgeMesh(pBEdgeMesh);
}

uint CMeshFactory::getNumOfBounaryEdgeMesh(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    return mpTMesh->getNumOfBoundaryEdgeMesh();
}

// "BoundaryNodeMesh"の節点(BoundarySBNode)生成
// ----
void CMeshFactory::GeneBoundaryNode(const uint& mgLevel, const uint& bnd_id, const uint& bndType,
                          const uint& mesh_id, const uint& node_id,
                          const uint& b_node_id, const uint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    CBoundaryNodeMesh* pBndNodeMesh= mpTMesh->getBndNodeMeshID(bnd_id);
    pBndNodeMesh->setBndType(bndType);
    
    // *DOF違いの同一BNodeがセットされるので,BNodeIDによる選別.
    // ----
    // BoundaryMeshに既に存在するか判定
    // ----
    uint crIndex= pBndNodeMesh->getNumOfBNode();
    uint crBNodeID;
    CBoundarySBNode *pCrBNode;
    uint ib; bool bfind(false);

    if(crIndex > 0){
        if(crIndex==1){
            pCrBNode= pBndNodeMesh->getBNodeIX(0);
            crBNodeID= pCrBNode->getID();

            if(b_node_id==crBNodeID) bfind= true;
            //cout << "crBNodeID= " << crBNodeID << ", b_node_id= " << b_node_id << endl;
        }
        for(ib=crIndex-1; ib > 0; ib--){
            pCrBNode= pBndNodeMesh->getBNodeIX(ib);
            crBNodeID= pCrBNode->getID();
            
            if(b_node_id==crBNodeID) bfind= true;
        };
    }
    
    
    if(bfind){
        // 既に生成済みのBNodeの場合
        pCrBNode= pBndNodeMesh->getBNodeID(b_node_id);
        pCrBNode->addDOF(dof);
        pCrBNode->setValue(dof, val);
    }else{
        // 新たにBNodeを生成する場合
        CNode *pNode= mpTMesh->getNode(node_id);
        
        CBoundarySBNode *pBNode = new CBoundarySBNode();

        pBNode->setNode(pNode);
        pBNode->setID(b_node_id);
        

        pBNode->addDOF(dof);
        pBNode->setValue(dof, val);
        pBndNodeMesh->addBNode(pBNode);
    }
}

// BoundaryFaceMeshの節点(BoundaryNode)生成
// ----
void CMeshFactory::GeneBoundaryFaceNode(const uint& mgLevel, const uint& bnd_id, const uint& bndType,
        const uint& mesh_id, const uint& node_id, const uint& b_node_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    // 境界面メッシュ
    CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshID(bnd_id);
    // ----------------------------------------------
    // 境界面メッシュへの属性設定は,GeneBoundaryFaceで実行
    // ----------------------------------------------


    CNode *pNode= mpTMesh->getNode(node_id);
    // 境界面メッシュのBNode
    CBoundaryNode *pBNode= new CBoundaryNode;// <<<<<<<<<<<<<< new
    pBNode->setNode(pNode);
    pBNode->setID(b_node_id);
    pBNode->setMGLevel(mgLevel);
    pBNode->resizeValue(mMGLevel+1);// 階層数+1  ,2010.07.28

    pBFaceMesh->addBNode(pBNode);//面メッシュへBNodeを追加
}
// BNode周辺のAggregate情報の初期化
//
void CMeshFactory::initFaceAggregate(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    // 境界面メッシュ
    CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshID(bnd_id);

    pBFaceMesh->resizeAggFace();
    pBFaceMesh->setupAggFace();
}
// BoundaryFaceMeshの面(BoundaryFace) 生成
// ----
void CMeshFactory::GeneBoundaryFace(const uint& mgLevel, const uint& bnd_id, const uint& bndType, const uint& elemType,
                            const uint& mesh_id, const uint& elem_id, const uint& face_id, vuint& vBNodeID,
                            const uint& b_face_id, const uint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    // 境界面メッシュ
    CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshID(bnd_id);
    
    pBFaceMesh->setBndType(bndType);
    pBFaceMesh->setID(bnd_id);

    //cout << "MeshFactory::GeneBoundaryFace A, pBFaceMesh_ID = " << pBFaceMesh->getID() << endl;

    // 境界面
    // *DOF違いの同一BFaceがセットされるので,BFaceIDによる選別.
    // ----
    // BoundaryMeshに既に存在するか判定
    // ----
    uint crIndex= pBFaceMesh->getNumOfBFace();
    uint crBFaceID;
    CBoundaryFace *pCrBFace;
    uint ib; bool bfind(false);
    if(crIndex > 0){
        if(crIndex==1){
            pCrBFace= pBFaceMesh->getBFaceIX(0);
            crBFaceID= pCrBFace->getID();

            if(b_face_id==crBFaceID) bfind= true;
        }
        for(ib=crIndex-1; ib > 0; ib--){
            pCrBFace= pBFaceMesh->getBFaceIX(ib);
            crBFaceID= pCrBFace->getID();

            if(b_face_id==crBFaceID) bfind= true;
        };
    }

    CBoundaryFace *pBFace;
    
    if(bfind){
        // 既存のFaceMeshにBFaceが存在する場合
        pBFace= pBFaceMesh->getBFaceID(b_face_id);
        pBFace->setBndValue(dof,val);

        //cout << "pBFace id= " << pBFace->getID() << endl;
        
    }else{
        // 新規にBFaceを生成する場合
        pBFace = new CBoundaryFace();// <<<<<<<<<<<<<< new
        
        pBFace->setElementID(elem_id);
        pBFace->setElementFaceID(face_id);
        pBFace->setID(b_face_id);
        pBFace->setBFaceShape(elemType);
        //pBFace->addDOF(dof);
        pBFace->setBndValue(dof,val);

        CElement *pElem= mpTMesh->getElement(elem_id);
        pBFace->setElement(pElem);
        
        //cout << "pBFace id= " << pBFace->getID() << endl;

        switch(elemType){
            case(ElementType::Quad):case(ElementType::Quad2):
                pBFace->resizeBNode(4);
                break;
            case(ElementType::Triangle):case(ElementType::Triangle2):
                pBFace->resizeBNode(3);
                break;
            default:
                //TODO:Logger
                break;
        }

        pBFaceMesh->addBFace(pBFace);
    }

    //新規にBFaceを生成した場合のみBNodeをセット
    if(!bfind){
        // 境界面節点
        CBoundaryNode *pBNode;
        uint numOfVert= vBNodeID.size();
        uint ivert;
        for(ivert=0; ivert < numOfVert; ivert++){

            pBNode= pBFaceMesh->getBNodeID(vBNodeID[ivert]);
            pBFace->setBNode(ivert, pBNode);

            // BoundaryMesh側でDOF管理に移行
            //    //DOFのセット
            //    // *重複する可能性があるので,現状のDOFをチェックの上でdof追加.
            //    uint numOfDOF= pBNode->getNumOfDOF();
            //    uint idof, crDOF;
            //    bool bCheck(false);
            //    for(idof=0; idof < numOfDOF; idof++){
            //        crDOF= pBNode->getDOF(idof);
            //        if(crDOF==dof) bCheck= true;
            //    };
            //    if(!bCheck) pBNode->addDOF(dof);
        };

        pBFace->calcArea();
    }
}

// BoundaryVolumeMeshの節点(BoundaryNode)生成
// ----
void CMeshFactory::GeneBoundaryVolumeNode(const uint& mgLevel, const uint& bnd_id, const uint& bndType,
        const uint& mesh_id, const uint& node_id, const uint& b_node_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    // 境界体積メッシュ
    CBoundaryVolumeMesh *pBVolumeMesh= mpTMesh->getBndVolumeMeshID(bnd_id);
    // ----------------------------------------------
    // 境界体積メッシュへの属性設定は,GeneBoundaryVolumeで実行
    // ----------------------------------------------

    CNode *pNode= mpTMesh->getNode(node_id);
    // 境界面メッシュのBNode
    CBoundaryNode *pBNode= new CBoundaryNode;// <<<<<<<<<<<<< new
    pBNode->setNode(pNode);
    pBNode->setID(b_node_id);
    pBNode->setMGLevel(mgLevel);
    pBNode->resizeValue(mMGLevel+1);// 階層数+1 ,2010.07.28

    pBVolumeMesh->addBNode(pBNode);//体積メッシュへBNodeを追加
}
// BNode周囲のAggregateVolume情報の初期化
//
void CMeshFactory::initVolumeAggregate(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    // 境界体積メッシュ
    CBoundaryVolumeMesh *pBVolumeMesh= mpTMesh->getBndVolumeMeshID(bnd_id);

    pBVolumeMesh->resizeAggVol();
    pBVolumeMesh->setupAggVol();
}
// BoundaryVolumeMeshの体積(BoundaryVolume) 生成
// ----
void CMeshFactory::GeneBoundaryVolume(const uint& mgLevel, const uint& bnd_id, const uint& bndType, const uint& elemType,
                            const uint& mesh_id, const uint& elem_id, vuint& vBNodeID,
                            const uint& b_vol_id, const uint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    // 境界体積メッシュ
    CBoundaryVolumeMesh *pBVolumeMesh= mpTMesh->getBndVolumeMeshID(bnd_id);
    pBVolumeMesh->setBndType(bndType);
    pBVolumeMesh->setID(bnd_id);

    // 境界体積
    // *DOF違いの同一BVolumeがセットされるので,BVolumeIDによる選別.
    // ----
    // BoundaryMeshに既に存在するか判定
    // ----
    uint crIndex= pBVolumeMesh->getNumOfVolume();
    uint crBVolumeID;
    CBoundaryVolume *pCrBVolume;
    uint ib; bool bfind(false);
    if(crIndex > 0){
        if(crIndex==1){
            pCrBVolume= pBVolumeMesh->getBVolumeIX(0);
            crBVolumeID= pCrBVolume->getID();

            if(b_vol_id==crBVolumeID) bfind= true;
        }
        for(ib=crIndex-1; ib > 0; ib--){
            pCrBVolume= pBVolumeMesh->getBVolumeIX(ib);
            crBVolumeID= pCrBVolume->getID();

            if(b_vol_id==crBVolumeID) bfind= true;
        };
    }

    CBoundaryVolume *pBVolume;

    if(bfind){
        // 既存のVolumeMeshにBVolumeが存在する場合
        pBVolume= pBVolumeMesh->getBVolumeID(b_vol_id);
        //pBVolume->addDOF(dof);
        pBVolume->setBndValue(dof,val);

    }else{
        // 新規にBVolumeを生成する場合
        // ----
        switch(elemType){
            case(ElementType::Hexa):case(ElementType::Hexa2):
                pBVolume = new CBoundaryHexa;// <<<<<<<<<<<<<< new
                break;
            case(ElementType::Tetra):case(ElementType::Tetra2):
                pBVolume = new CBoundaryTetra;// <<<<<<<<<<<<< new
                break;
            case(ElementType::Prism):case(ElementType::Prism2):
                pBVolume = new CBoundaryPrism;// <<<<<<<<<<<<< new
                break;
            default:
                mpLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CMeshFactory::GeneBoundaryVolume");
                break;
        }

        pBVolume->setElementID(elem_id);
        pBVolume->setID(b_vol_id);
        pBVolume->setBndValue(dof,val);

        CElement *pElem= mpTMesh->getElement(elem_id);
        pBVolume->setElement(pElem);

        switch(elemType){
            case(ElementType::Hexa):case(ElementType::Hexa2):
                pBVolume->resizeBNode(8);
                break;
            case(ElementType::Tetra):case(ElementType::Tetra2):
                pBVolume->resizeBNode(4);
                break;
            case(ElementType::Prism):case(ElementType::Prism2):
                pBVolume->resizeBNode(6);
                break;
            default:
                // TODO:Logger
                break;
        }

        pBVolumeMesh->addBVolume(pBVolume);
    }

    //新規のBVolumeを生成した場合のみBNodeをセット
    if(!bfind){
        // 境界体積節点
        CBoundaryNode *pBNode;
        uint numOfVert= vBNodeID.size();
        uint ivert;
        for(ivert=0; ivert < numOfVert; ivert++){
            pBNode= pBVolumeMesh->getBNodeID(vBNodeID[ivert]);
            pBVolume->setBNode(ivert, pBNode);

            // BoundaryMeshでDOF管理に移行
            // ---------
            //    //DOFのセット
            //    // *重複する可能性があるので,現状のDOFをチェックの上でdof追加.
            //    uint numOfDOF= pBNode->getNumOfDOF();
            //    uint idof, crDOF;
            //    bool bCheck(false);
            //    for(idof=0; idof < numOfDOF; idof++){
            //        crDOF= pBNode->getDOF(idof);
            //        if(crDOF==dof) bCheck= true;
            //    };
            //    if(!bCheck) pBNode->addDOF(dof);
        };

        pBVolume->calcVolume();
    }
}

// BoundaryEdgeMeshの節点(BoundaryNode)生成
// ----
void CMeshFactory::GeneBoundaryEdgeNode(const uint& mgLevel, const uint& bnd_id, const uint& bndType,
        const uint& mesh_id, const uint& node_id, const uint& b_node_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    // 境界辺メッシュ
    CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshID(bnd_id);
    // ----------------------------------------------
    // 境界辺メッシュへの属性設定は,GeneBoundaryEdgeで実行
    // ----------------------------------------------

    CNode *pNode= mpTMesh->getNode(node_id);
    // 境界面メッシュのBNode
    CBoundaryNode *pBNode= new CBoundaryNode;// <<<<<<<<<<<<< new


    pBNode->setNode(pNode);
    pBNode->setID(b_node_id);
    pBNode->setMGLevel(mgLevel);
    pBNode->resizeValue(mMGLevel+1);// 階層数+1 ,2010.07.28


    ////debug
    //cout << "MeshFactory::GeneBoundaryEdgeNode, BNodeID= " << pBNode->getID() << ", NodeID= " << pNode->getID() << endl;

    pBEdgeMesh->addBNode(pBNode);//辺メッシュへBNodeを追加
}

// BNode周囲のAggregateEdge情報の初期化
//
void CMeshFactory::initEdgeAggregate(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    // 境界辺メッシュ
    CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshID(bnd_id);

    pBEdgeMesh->resizeAggEdge();
    pBEdgeMesh->setupAggEdge();
}

// BoundaryEdgeMeshの辺(BoundaryEdge) 生成
// ----
void CMeshFactory::GeneBoundaryEdge(const uint& mgLevel, const uint& bnd_id, const uint& bndType, const uint& elemType,
        const uint& mesh_id, const uint& elem_id, const uint& edge_id, vuint& vBNodeID,
        const uint& b_edge_id, const uint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    // 境界辺メッシュ
    CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshID(bnd_id);
    pBEdgeMesh->setBndType(bndType);
    pBEdgeMesh->setID(bnd_id);

    // 境界辺
    // *DOF違いの同一BEdgeがセットされるので,BEdgeIDによる選別.
    // ----
    // カレントBEdgeIDを取得
    // ----
    uint crIndex= pBEdgeMesh->getNumOfEdge();
    uint crBEdgeID;
    CBoundaryEdge *pCrBEdge;
    uint ib;
    bool bfind(false);//新規生成か否かのフラグ
    if(crIndex > 0){
        if(crIndex==1){
            pCrBEdge= pBEdgeMesh->getBEdgeIX(0);
            crBEdgeID= pCrBEdge->getID();

            if(b_edge_id==crBEdgeID) bfind= true;
        }
        for(ib=crIndex-1; ib > 0; ib--){
            pCrBEdge= pBEdgeMesh->getBEdgeIX(ib);
            crBEdgeID= pCrBEdge->getID();

            if(b_edge_id==crBEdgeID) bfind= true;
        };
    }

    CBoundaryEdge *pBEdge;
    
    if(bfind){
        // 既存のEdgeMeshにBEdgeが存在する場合
        pBEdge= pBEdgeMesh->getBEdgeID(b_edge_id);
        //pBEdge->addDOF(dof);
        pBEdge->setBndValue(dof,val);

    }else{
        //新規にBEdgeを生成する場合
        pBEdge = new CBoundaryEdge();// <<<<<<<<<<<<<< new

        pBEdge->setElementID(elem_id);
        pBEdge->setElementEdgeID(edge_id);
        pBEdge->setID(b_edge_id);
        pBEdge->setBndValue(dof,val);
        pBEdge->setBEdgeShape(elemType);

        CElement *pElem= mpTMesh->getElement(elem_id);
        pBEdge->setElement(pElem);

        switch(elemType){
            case(ElementType::Beam):case(ElementType::Beam2):
                pBEdge->resizeBNode(2);
                break;
            default:
                //TODO:Logger
                break;
        }

        pBEdgeMesh->addBEdge(pBEdge);
    }

    if(!bfind){
        // 境界辺節点
        CBoundaryNode *pBNode;
        uint numOfVert= vBNodeID.size();
        uint ivert;
        for(ivert=0; ivert < numOfVert; ivert++){
            pBNode= pBEdgeMesh->getBNodeID(vBNodeID[ivert]);
            pBEdge->setBNode(ivert, pBNode);
            
            // BoundaryMeshでDOF管理に移行.
            // ---
            //    //DOFのセット
            //    // *重複する可能性があるので,現状のDOFをチェックの上でdof追加.
            //    uint numOfDOF= pBNode->getNumOfDOF();
            //    uint idof, crDOF;
            //    bool bCheck(false);
            //    for(idof=0; idof < numOfDOF; idof++){
            //        crDOF= pBNode->getDOF(idof);
            //        if(crDOF==dof) bCheck= true;
            //    };
            //    if(!bCheck) pBNode->addDOF(dof);
        };

        pBEdge->calcLength();
    }
}

// ----
// 境界条件階層化
// ----
void CMeshFactory::refineBoundary()
{
    uint numOfMesh;
    uint iLevel, iMesh;
    for(iLevel=0; iLevel < mMGLevel; iLevel++){

        mpTAssyModel = mpGMGModel->getAssyModel(iLevel);
        CAssyModel *pProgAssy= mpGMGModel->getAssyModel(iLevel+1);

        numOfMesh= mpTAssyModel->getNumOfMesh();
        for(iMesh=0; iMesh < numOfMesh; iMesh++){
            
            mpTMesh = mpTAssyModel->getMesh(iMesh);
            CMesh *pProgMesh= pProgAssy->getMesh(iMesh);

            // BoundaryFaceMesh
            // ----
            uint numOfBFaceMesh= mpTMesh->getNumOfBoundaryFaceMesh();
            pProgMesh->reserveBndFaceMesh(numOfBFaceMesh);
            uint iBFaceMesh;
            for(iBFaceMesh=0; iBFaceMesh < numOfBFaceMesh; iBFaceMesh++){
                
                CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshIX(iBFaceMesh);
                
                CBoundaryFaceMesh *pProgBFaceMesh= new CBoundaryFaceMesh;// <<<<<<<<<<<<<<<<< new
                
                pProgMesh->setBndFaceMesh(pProgBFaceMesh);

                //progBFaceMeshへDOFをセット
                pProgBFaceMesh->resizeDOF(pBFaceMesh->getNumOfDOF());
                uint idof, dof;
                for(idof=0; idof < pBFaceMesh->getNumOfDOF(); idof++){
                    dof = pBFaceMesh->getDOF(idof);
                    pProgBFaceMesh->setDOF(idof, dof);
                }

                pBFaceMesh->GeneEdgeBNode();
                pBFaceMesh->GeneFaceBNode();
                
                pBFaceMesh->refine(pProgBFaceMesh);// Refine
                pProgBFaceMesh->resizeAggFace();
                pProgBFaceMesh->setupAggFace();
                
                uint nBndType= pBFaceMesh->getBndType();
                pProgBFaceMesh->setBndType(nBndType);// Dirichlet.or.Neumann

                uint nID= pBFaceMesh->getID();
                pProgBFaceMesh->setID(nID);

                pProgBFaceMesh->setMGLevel(iLevel+1);
                pProgBFaceMesh->setMaxMGLevel(mMGLevel);
                
                //BNodeへの境界値の再配分
                //--
                if(BoundaryType::Dirichlet==pBFaceMesh->getBndType()){
                    pBFaceMesh->distValueBNode();//Dirichlet:親の方でBNodeへの配分を決める
                }
                if(BoundaryType::Neumann==pBFaceMesh->getBndType()){
                    if(iLevel==0) pBFaceMesh->distValueBNode();//Level=0の場合は、親側でNeumannも分配
                    pProgBFaceMesh->distValueBNode();//Neumann:子供の方でBNOdeへの配分を決める
                }
            };

            // BoundaryEdgeMesh
            // ----
            uint numOfEdgeMesh= mpTMesh->getNumOfBoundaryEdgeMesh();
            pProgMesh->reserveBndEdgeMesh(numOfEdgeMesh);
            uint iBEdgeMesh;
            for(iBEdgeMesh=0; iBEdgeMesh < numOfEdgeMesh; iBEdgeMesh++){

                CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshIX(iBEdgeMesh);
                
                CBoundaryEdgeMesh *pProgBEdgeMesh= new CBoundaryEdgeMesh;// <<<<<<<<<<<<<< new
                
                pProgMesh->setBndEdgeMesh(pProgBEdgeMesh);

                //ProgBEdgeMeshにDOFをセット
                pProgBEdgeMesh->resizeDOF(pBEdgeMesh->getNumOfDOF());
                uint idof, dof;
                for(idof=0; idof < pBEdgeMesh->getNumOfDOF(); idof++){
                    dof = pBEdgeMesh->getDOF(idof);
                    pProgBEdgeMesh->setDOF(idof, dof);
                }

                pBEdgeMesh->GeneEdgeBNode();
                pBEdgeMesh->refine(pProgBEdgeMesh);// Refine

                pProgBEdgeMesh->resizeAggEdge();
                pProgBEdgeMesh->setupAggEdge();

                uint nBndType= pBEdgeMesh->getBndType();
                pProgBEdgeMesh->setBndType(nBndType);// Dirichlet.or.Neumann

                uint nID= pBEdgeMesh->getID();
                pProgBEdgeMesh->setID(nID);

                pProgBEdgeMesh->setMGLevel(iLevel+1);
                pProgBEdgeMesh->setMaxMGLevel(mMGLevel);

                //BNodeへの境界値の再配分
                //--
                if(BoundaryType::Dirichlet==pBEdgeMesh->getBndType()){
                    pBEdgeMesh->distValueBNode();//Dirichlet:親の方でBNodeへの配分を決める
                }
                if(BoundaryType::Neumann==pBEdgeMesh->getBndType()){
                    if(iLevel==0) pBEdgeMesh->distValueBNode();//Level=0の場合は、親側でNeumannも分配
                    pProgBEdgeMesh->distValueBNode();//Neumann:子供の方でBNOdeへの配分を決める
                }
            };

            // BoundaryVolumeMesh
            // ----
            uint numOfVolMesh= mpTMesh->getNumOfBoundaryVolumeMesh();
            pProgMesh->reserveBndVolumeMesh(numOfVolMesh);
            uint iBVolMesh;
            for(iBVolMesh=0; iBVolMesh < numOfVolMesh; iBVolMesh++){

                CBoundaryVolumeMesh *pBVolMesh= mpTMesh->getBndVolumeMeshIX(iBVolMesh);
                
                CBoundaryVolumeMesh *pProgBVolMesh= new CBoundaryVolumeMesh;// <<<<<<<<<<<<<<<<<< new

                pProgMesh->setBndVolumeMesh(pProgBVolMesh);

                //pProgBVolMeshにDOFをセット
                pProgBVolMesh->resizeDOF(pBVolMesh->getNumOfDOF());
                uint idof, dof;
                for(idof=0; idof < pBVolMesh->getNumOfDOF(); idof++){
                    dof = pBVolMesh->getDOF(idof);
                    pProgBVolMesh->setDOF(idof, dof);
                }

                pBVolMesh->GeneEdgeBNode();
                pBVolMesh->GeneFaceBNode();
                pBVolMesh->GeneVolBNode();
                
                pBVolMesh->refine(pProgBVolMesh);// Refine
                pProgBVolMesh->resizeAggVol();
                pProgBVolMesh->setupAggVol();

                uint nBndType= pBVolMesh->getBndType();
                pProgBVolMesh->setBndType(nBndType);// Dirichlet .or. Neumann

                uint nID= pBVolMesh->getID();
                pProgBVolMesh->setID(nID);

                pProgBVolMesh->setMGLevel(iLevel+1);
                pProgBVolMesh->setMaxMGLevel(mMGLevel);
                
                //BNodeへの境界値の再配分
                //--
                if(BoundaryType::Dirichlet==pBVolMesh->getBndType()){
                    pBVolMesh->distValueBNode();//Dirichlet:親の方でBNodeへの配分を決める
                }
                if(BoundaryType::Neumann==pBVolMesh->getBndType()){
                    if(iLevel==0) pBVolMesh->distValueBNode();//Level=0の場合は、親側でNeumannも分配
                    pProgBVolMesh->distValueBNode();//Neumann:子供の方でBNOdeへの配分を決める
                }
            };

            // BoundaryNodeMesh : BNodeMeshGrp* を各階層にセット.
            // ----
            CBNodeMeshGrp* pBNodeMeshGrp= mpTMesh->getBNodeMeshGrp();
            pProgMesh->setBNodeMeshGrp(pBNodeMeshGrp);
        };
    };
}



//----
// Bucket in Mesh
//----

// Node_Index set to Bucket(in Mesh):: All in One Method
//
void CMeshFactory::setupBucketNode(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->setupBucketNode();// Bucket一括処理

}
// initialize Bucket for Node
void CMeshFactory::initBucketNode(const uint& mgLevel, const uint& mesh_id, const uint& maxID, const uint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->initBucketNode(maxID, minID);// Bucket領域確保
}
// set (ID & Index) to Bucket for Node
void CMeshFactory::setIDBucketNode(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& index)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->setupBucketNodeIndex(id, index);

}

// Element_Index set to Bucket (in Mesh):: All in-One Method
//
void CMeshFactory::setupBucketElement(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->setupBucketElement();
}
// initialize Bucket for Element
void CMeshFactory::initBucketElement(const uint& mgLevel, const uint& mesh_id, const uint& maxID, const uint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    ;
    mpTMesh->initBucketElement(maxID, minID);
}
// set (ID & Index) to Bucket for Element
void CMeshFactory::setIDBucketElement(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& index)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->setupBucketElementIndex(id, index);
}

//void CMeshFactory::setMaxMinID_TMeshN(const uint &maxID_Node, const uint &minID_Node)
//{
//    mpTBucket->resizeBucketNode(maxID_Node, minID_Node);
//}
//
//// Element index set for mpTBucket
////
//void CMeshFactory::setMaxMinID_TMeshE(const uint &maxID_Elem, const uint &minID_Elem)
//{
//    mpTBucket->resizeBucketElement(maxID_Elem, minID_Elem);
//}


// 材質データの配列の確保
// --
void CMeshFactory::reserveMaterial(const uint& res_size)
{
     mpGMGModel->reserveMaterial(res_size);
}

// 材質データの生成
// --
void CMeshFactory::GeneMaterial(const uint& mesh_id, const uint& material_id, string& name, vuint& vType, vdouble& vValue)
{
    CMaterial *pMaterial = new CMaterial;

    pMaterial->setID(material_id);
    pMaterial->setName(name);
    pMaterial->setMeshID(mesh_id);

    for(uint i=0; i< vType.size(); i++) pMaterial->setValue(vType[i],vValue[i]);

    
    mpGMGModel->setMaterial(pMaterial);
}




// 通信領域(CommMesh)の配列確保
//
void CMeshFactory::reserveCommMesh(const uint& mgLevel, const uint& mesh_id, const uint& res_size)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->reserveCommMesh(res_size);
}

// 通信領域(CommMesh)の生成
//
void CMeshFactory::GeneCommMesh(const uint& mgLevel, const uint& mesh_id, const uint& comID, const uint& myRank, const uint& nTransmitRank)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    CIndexBucket *pBucket= mpTMesh->getBucket();

    mpTCommMesh = new CCommMesh(pBucket);
    mpTCommMesh->setCommID(comID);
    mpTCommMesh->setRankID(myRank);
    mpTCommMesh->setTransmitRankID(nTransmitRank);

    
    mpTMesh->setCommMesh(mpTCommMesh);
}

// 通信ノード(CommMesh内に,MeshのNodeポインターの配列を確保)
//  => CommNodeというものは,存在しない.
//
void CMeshFactory::reserveCommNode(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, const uint& res_size)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);

    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    
    mpTCommMesh->reserveNode(res_size);   //CommMeshノードの配列予約
    mpTCommMesh->resizeNodeRank(res_size);//NodeRank一覧の確保

}

// MeshからNodeを取得してCommMeshにセット
//  => CommNodeというものは,存在しない.
//
void CMeshFactory::GeneCommNode(const uint& mgLevel, const uint& commNodeID,
                                     const uint& mesh_id, const uint& commesh_id, const uint& nodeID, const uint& rank)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);

    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);

    // MeshからNodeを取得して,CommMeshにセット
    // --
    CNode* pNode= mpTMesh->getNode(nodeID);
    mpTCommMesh->setNode(pNode);// 順番にpush_backしている.
    mpTCommMesh->setNodeRank(commNodeID, rank);//Nodeランク一覧 配列

    // Send,Recvノードのセット
    if(mpTCommMesh->getRankID()==rank) mpTCommMesh->setSendNode(pNode, commNodeID);
    if(mpTCommMesh->getTransmitRankID()==rank) mpTCommMesh->setRecvNode(pNode, commNodeID);
}

// 通信要素(CommElement)
//
void CMeshFactory::reserveCommElement(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, const uint& res_size)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);

    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    mpTCommMesh->reserveCommElementAll(res_size);
}

// CommElementを生成して,MeshからElementを取得してセット
//
void CMeshFactory::GeneCommElement(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, 
                                   const uint& nType, const uint& elemID, vuint& vCommNodeID)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);

    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    
    uint numOfVert= vCommNodeID.size();
    vuint vNodeRank;  vNodeRank.reserve(numOfVert);
    
    uint ivert, rank;
    uint commNodeID;
    for(ivert=0; ivert< numOfVert; ivert++){
        commNodeID= vCommNodeID[ivert];
        rank= mpTCommMesh->getNodeRank(commNodeID);
        
        vNodeRank.push_back(rank);
    };
    
    CCommElement *pCommElem;
    switch(nType){
        case(ElementType::Hexa):
            pCommElem = new CCommHexa;
            break;
        case(ElementType::Tetra):
            pCommElem = new CCommTetra;
            break;
        case(ElementType::Prism):
            pCommElem = new CCommPrism;
            break;
//        case(ElementType::Pyramid):
//            pCommElem = new CCommPyramid;
//            break;
        case(ElementType::Quad):
            pCommElem = new CCommQuad;
            break;
        case(ElementType::Triangle):
            pCommElem = new CCommTriangle;
            break;
        case(ElementType::Beam):
            pCommElem = new CCommBeam;
            break;
    }
    
    CElement* pElem= mpTMesh->getElement(elemID);

    //cout << "GeneCommElement  pElem ID==" << pElem->getID() << endl;

    pCommElem->setElement(pElem);
    
    // Nodeランクのセット
    // --
    for(ivert=0; ivert< numOfVert; ivert++){
        rank = vNodeRank[ivert];
        pCommElem->setNodeRank(ivert, rank);
    };
    
    mpTCommMesh->setCommElementAll(pCommElem);
    
//
//    prolongation時の処理 Send,Recvの割り振りがファイル入力時に行われていないため参考.
//
//    //1.CommMesh内でのCommElemのIndex番号の割り振り && CommMesh内のCommElementを通信するか否かのCommElementに選別
//    pProgCommMesh->AllocateCommElement();
//
//    //2.CommMesh内でのCommElemの隣接情報
//    //3.CommMesh内でのNodeのIndex番号の割り振り,CommMeshのmvNode,mvSendNode,mvRecvNodeの取得
//    pProgCommMesh->setupAggCommElement(pProgMesh->getElements());
//    pProgCommMesh->sortCommNodeIndex();// CommMesh内でのNode Index番号生成, Send,Recvノードの選別, DNode,DElementの選別ソート
//                                       // mvNodeのセットアップもsortCommNodeIndexから,CommElementのsetCommNodeIndex内でmvNodeにセットしている.
//    //4.mapデータのセットアップ
//    pProgCommMesh->setupMapID2CommID();
}

//// 通信領域のprolongation
//// ○ refineMesh();後に呼ばれるルーチン
//// --
//void CMeshFactory::refineCommMesh()
//{
//    //-- Mesh
//    CAssyModel *pAssy, *pProgAssy;
//    CMesh      *pMesh, *pProgMesh;
//    //-- 通信Mesh
//    CCommMesh  *pCommMesh, *pProgCommMesh;
//    //-- 通信要素(CommElem)
//    CCommElement      *pCommElem;//元になったCommElem(親のCommElem)
//    vector<CCommElement*> vProgCommElem;//生成されるprogCommElem達
//    CCommElement         *pProgCommElem;//progCommElem <= prolongation CommElem
//
//
//    uint numOfMesh, numOfCommMesh, numOfCommElemAll, numOfCommNode;
//    uint ilevel,imesh,icommesh,icomelem,iprocom,ivert;
//    // ---
//    // 階層Level ループ
//    // ---
//    for(ilevel=0; ilevel< mMGLevel; ilevel++){
//
//        //debug
//        cout << "ilevel => " << ilevel << endl;
//
//        pAssy= mpGMGModel->getAssyModel(ilevel);
//        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);
//
//        numOfMesh= pAssy->getNumOfMesh();
//
//        ////debug
//        //cout << "numOfMesh =>" << numOfMesh << endl;
//
//        // ---
//        // Mesh(パーツ) ループ in AssyMode
//        // ---
//        for(imesh=0; imesh< numOfMesh; imesh++){
//            pMesh= pAssy->getMesh(imesh);
//            pProgMesh= pProgAssy->getMesh(imesh);
//
//            numOfCommMesh= pMesh->getNumOfCommMesh();
//            ////debug
//            //cout << "numOfCommMesh => " << numOfCommMesh << endl;
//
//            // --
//            // CommMesh(通信領域) ループ in Mesh
//            // --
//            for(icommesh=0; icommesh< numOfCommMesh; icommesh++){
//                pCommMesh= pMesh->getCommMesh(icommesh);
//
//                // "new CommMesh" に下段階層のプロパティ・セット
//                // --
//                pProgCommMesh= new CCommMesh;// <<<<<<<<<<<-- prolongation CommMesh
//                pProgCommMesh->setCommID( pCommMesh->getCommID());
//                pProgCommMesh->setRankID( pCommMesh->getRankID());
//                pProgCommMesh->setTransmitRankID( pCommMesh->getTransmitRankID());
//
//                pProgMesh->setCommMesh(pProgCommMesh);// プロパティを全てセットしてからMeshにCommMeshをセットすること.
//
//
//                numOfCommElemAll= pCommMesh->getNumOfCommElementAll();
//                pProgCommMesh->reserveCommElementAll(numOfCommElemAll*8);// <<<<<<-- CommElemAllリザーブ
//
//                numOfCommNode= pCommMesh->getNumOfNode();
//                pProgCommMesh->reserveNode(numOfCommNode*8);// <<<<<<<<<<<<-- CommMeshノード リザーブ
//
//                //debug
//                string sOutputStr =  boost::lexical_cast<string>(pCommMesh->getCommID()) + ","
//                                    +boost::lexical_cast<string>(pCommMesh->getRankID()) + ","
//                                    +boost::lexical_cast<string>(pCommMesh->getTransmitRankID()) + ","
//                                    +boost::lexical_cast<string>(numOfCommElemAll);
//                //debug
//                mpLogger->Info(Utility::LoggerMode::Debug,"setup to progCommMesh:commID,rank,transmit_rank,CommElemAll数 => ",sOutputStr);
//
//
//                //pCommMeshからCommElemを取得し,progCommElemを生成 => progCommMeshにセット
//                // 計算rank(計算領域番号)は,CommMeshが所有.
//                for(icomelem=0; icomelem< numOfCommElemAll; icomelem++){
//
//                    pCommElem= pCommMesh->getCommElementAll(icomelem);
//                    pCommElem->setupProgNodeRank(ilevel+1);//辺,面,体積中心にRank設定. <<<<<<-- progCommElemのNodeランク
//
//                    vProgCommElem.clear();
//
//                    //debug
//                    mpLogger->Info(Utility::LoggerMode::Debug,"commelem index=>", icomelem);
//
//                    // 形状別のprogCommElem
//                    switch(pCommElem->getShapeType()){
//                        case(ElementType::Hexa):
//                            vProgCommElem.reserve(8);
//                            for(ivert=0; ivert< 8; ivert++){
//                                pProgCommElem= new CCommHexa;// <<<<<<<<<-- prolongation CommElement
//                                vProgCommElem.push_back(pProgCommElem);
//                            };
//                            dividCommElem(pCommElem, vProgCommElem);
//
//                            //debug
//                            mpLogger->Info(Utility::LoggerMode::Debug,"CommElem(Hexa)の分割");
//
//                            break;
//                        case(ElementType::Tetra):
//                            vProgCommElem.reserve(4);
//                            for(ivert=0; ivert< 4; ivert++){
//                                pProgCommElem= new CCommHexa;// <<<<<<<<<-- prolongation CommElement
//                                vProgCommElem.push_back(pProgCommElem);
//                            };
//                            dividCommElem(pCommElem, vProgCommElem);
//
//                            break;
//                        case(ElementType::Prism):
//                            vProgCommElem.reserve(6);
//                            for(ivert=0; ivert< 6; ivert++){
//                                pProgCommElem= new CCommHexa;// <<<<<<<<<-- prolongation CommElement
//                                vProgCommElem.push_back(pProgCommElem);
//                            };
//                            dividCommElem(pCommElem, vProgCommElem);
//
//                            break;
//                        case(ElementType::Pyramid):
//                            // Pyramid => Hexa縮退 ?
//                            vProgCommElem.reserve(8);
//                            for(ivert=0; ivert< 8; ivert++){
//                                pProgCommElem= new CCommHexa;// <<<<<<<<<-- prolongation CommElement
//                                vProgCommElem.push_back(pProgCommElem);
//                            };
//                            dividCommElem(pCommElem, vProgCommElem);
//
//                            break;
//                        case(ElementType::Quad):
//                            vProgCommElem.reserve(4);
//                            for(ivert=0; ivert< 4; ivert++){
//                                pProgCommElem= new CCommQuad;// <<<<<<<<<-- prolongation CommElement
//                                vProgCommElem.push_back(pProgCommElem);
//                            };
//                            dividCommElem(pCommElem, vProgCommElem);
//
//                            break;
//                        case(ElementType::Triangle):
//                            vProgCommElem.reserve(3);
//                            for(ivert=0; ivert< 3; ivert++){
//                                pProgCommElem= new CCommQuad;// <<<<<<<<<-- prolongation CommElement
//                                vProgCommElem.push_back(pProgCommElem);
//                            };
//                            dividCommElem(pCommElem, vProgCommElem);
//
//                            break;
//                        case(ElementType::Beam):
//                            vProgCommElem.reserve(2);
//                            for(ivert=0; ivert< 2; ivert++){
//                                pProgCommElem= new CCommBeam;// <<<<<<<<<-- prolongation CommElement
//                                vProgCommElem.push_back(pProgCommElem);
//                            };
//                            dividCommElem(pCommElem, vProgCommElem);
//
//                            break;
//                        default:
//                            mpLogger->Info(Utility::LoggerMode::Error, "refineCommMesh内の,CommElement ShapeType Error @MeshFactory::refineCommMesh");
//                            break;
//                    }
//                    // progCommMeshへprogCommElmeAllのセット(全てのCommElement)
//                    // --
//                    for(iprocom=0; iprocom< vProgCommElem.size(); iprocom++){
//                        pProgCommMesh->setCommElementAll(vProgCommElem[iprocom]);
//                    };
//
//                    //debug
//                    mpLogger->Info(Utility::LoggerMode::Debug, "refineCommMesh内の,pProgCommMesh->setCommElem");
//
//                };//CommElem ループ
//
//                // 1.CommMesh内でのCommElemのIndex番号の割り振り && CommMesh内のCommElementを通信するか否かのCommElementに選別
//                pProgCommMesh->AllocateCommElement();
//
//                //debug
//                cout << "pProgCommMesh->AllocateCommElement" << endl;
//
//                // 2.CommMesh内でのCommElemの隣接情報
//                // 3.CommMesh内でのNodeのIndex番号の割り振り,CommMeshのmvNode,mvSendNode,mvRecvNodeの取得
//                //
//                pProgCommMesh->setupAggCommElement(pProgMesh->getElements());
//                //debug
//                cout << "pProgCommMesh->setupAggCommElement" << endl;
//
//                pProgCommMesh->sortCommNodeIndex();// CommMesh内でのNode Index番号生成, Send,Recvノードの選別, DNode,DElementの選別ソート
//                                                   // mvNodeのセットアップもsortCommNodeIndexから,CommElementのsetCommNodeIndex内でmvNodeにセットしている.
//                //debug
//                cout << "pProgCommMesh->sortCommNodeIndex" << endl;
//
//                // 4.mapデータのセットアップ
//                pProgCommMesh->setupMapID2CommID();
//
//                //debug
//                mpLogger->Info(Utility::LoggerMode::Debug, "refineCommMesh内の,pProgCommMesh->Allocate..setupMapID2CommID");
//
//            };//CommMesh ループ
//
//            // Mesh のNode,Elementの計算領域整理
//            // --
//            pProgMesh->sortMesh();//MeshのmvNode,mvElementから計算に使用しないNode(DNode),Element(DElement)を移動
//
//
//        };//Meshループ
//    };//Levelループ
//}

// prolongation CommElementの生成
// --
void CMeshFactory::GeneProgCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem)
{
    CCommElement *pProgCommElem;
    uint ivert;
    
    // 形状別のprogCommElem
    switch(pCommElem->getShapeType()){
        case(ElementType::Hexa):
            vProgCommElem.reserve(8);
            for(ivert=0; ivert< 8; ivert++){
                pProgCommElem= new CCommHexa;// <<<<<<<<<-- prolongation CommElement
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);

            //debug
            mpLogger->Info(Utility::LoggerMode::Debug,"CommElem(Hexa)の分割");

            break;
        case(ElementType::Tetra):
            vProgCommElem.reserve(4);
            for(ivert=0; ivert< 4; ivert++){
                pProgCommElem= new CCommHexa;// <<<<<<<<<-- prolongation CommElement
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);

            break;
        case(ElementType::Prism):
            vProgCommElem.reserve(6);
            for(ivert=0; ivert< 6; ivert++){
                pProgCommElem= new CCommHexa;// <<<<<<<<<-- prolongation CommElement
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);

            break;
//        case(ElementType::Pyramid):
//            // Pyramid => Hexa縮退 ?
//            vProgCommElem.reserve(8);
//            for(ivert=0; ivert< 8; ivert++){
//                pProgCommElem= new CCommHexa;// <<<<<<<<<-- prolongation CommElement
//                vProgCommElem.push_back(pProgCommElem);
//            };
//            dividCommElem(pCommElem, vProgCommElem);
//
//            break;
        case(ElementType::Quad):
            vProgCommElem.reserve(4);
            for(ivert=0; ivert< 4; ivert++){
                pProgCommElem= new CCommQuad;// <<<<<<<<<-- prolongation CommElement
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);

            break;
        case(ElementType::Triangle):
            vProgCommElem.reserve(3);
            for(ivert=0; ivert< 3; ivert++){
                pProgCommElem= new CCommQuad;// <<<<<<<<<-- prolongation CommElement
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);

            break;
        case(ElementType::Beam):
            vProgCommElem.reserve(2);
            for(ivert=0; ivert< 2; ivert++){
                pProgCommElem= new CCommBeam;// <<<<<<<<<-- prolongation CommElement
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);

            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error, "ShapeType Error @MeshFactory::GeneProgCommElem");
            break;
    }
}

// 1.再分割CommElement(progCommElem)へ,再分割Elementをセット
// 2.再分割CommElement(progCommElem)の頂点へ"rank"をセット
// --
void CMeshFactory::dividCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem)
{
    CProgElementTree *pProgTree = CProgElementTree::Instance();

    CElement* pElem= pCommElem->getElement();
    CElement* pProgElem;
    CCommElement* pProgCommElem;
    uint nRank;
    uint ivert, progvert;
    uint iedge, iface;
    uint invalid= pProgTree->getInvalidNum();

    ////debug
    //mpLogger->Info(Utility::LoggerMode::Debug,"Factory::dividCommElemの invalid => ",invalid);

    uint numOfVert, numOfEdge, numOfFace;
    numOfVert= pElem->getNumOfNode(); numOfEdge= pElem->getNumOfEdge(); numOfFace= pElem->getNumOfFace();

    // 親CommElemの頂点ループ(子CommElemのアドレス)
    for(ivert=0; ivert< numOfVert; ivert++){
        pProgElem= pElem->getProgElem(ivert);

        pProgCommElem= vProgCommElem[ivert];
        pProgCommElem->setElement(pProgElem);// progCommElemへ,子要素(progElment)のセット

        // progCommElemの各頂点へ,rank(DomID)のセット
        // --

        // 親のvertexのrankを子のvertexへセット
        progvert= pProgTree->getVertProgVert(ivert, pCommElem->getShapeType());
        nRank= pCommElem->getNodeRank(ivert);
        pProgCommElem->setNodeRank(progvert,nRank);
        ////debug
        //mpLogger->Info(Utility::LoggerMode::Debug,"Factory::dividCommElemのVertexRank=>setNodeRank",(uint)nRank);


        // edgeノードに対応するprogCommElemのrank
        for(iedge=0; iedge< numOfEdge; iedge++){
            
            progvert= pProgTree->getEdgeProgVert(iedge, ivert, pCommElem->getShapeType());

            if(progvert != invalid){
                nRank= pCommElem->getEdgeRank(iedge);
                pProgCommElem->setNodeRank(progvert, nRank);

                ////debug
                //mpLogger->Info(Utility::LoggerMode::Debug,"Factory::dividCommElemのEdgeRank=>setNodeRank",(uint)nRank);
            }
        };
        
        
        // faceノードに対応するprogCommElemのrank
        for(iface=0; iface< numOfFace; iface++){
            
            progvert= pProgTree->getFaceProgVert(iface, ivert, pCommElem->getShapeType());

            if(progvert != invalid){
                nRank= pCommElem->getFaceRank(iface);
                pProgCommElem->setNodeRank(progvert, nRank);

                ////debug
                //mpLogger->Info(Utility::LoggerMode::Debug,"Factory::dividCommElemのFaceRank=>setNodeRank",(uint)nRank);
            }
        };
        // volumeに対応するprogCommElemのrank
        progvert= pProgTree->getVolProgVert(ivert, pCommElem->getShapeType());
        nRank= pCommElem->getVolRank();
        pProgCommElem->setNodeRank(progvert, nRank);
        ////debug
        //mpLogger->Info(Utility::LoggerMode::Debug,"Factory::dividCommElemのVolRank=>setNodeRank",(uint)nRank);
    };
}




// [コンタクトメッシュ数ぶん呼び出される]
//  -------------------------------
// ファイル入力時にレベル0のContactMesh要素へのマーキング
//
//  -> 全レベルのContactMeshは,ここで予め生成しておく.
//  -> 自身と同じランクに所属する接合面は, setupContactMesh, setupSkin で生成される.
// --
// コンタクトメッシュを全階層に生成
// --
void CMeshFactory::GeneContactMesh(const uint& contactID, const uint& myRank, const uint& transRank, const uint& nProp)
{
    uint ilevel;
    for(ilevel=0; ilevel< mMGLevel+1; ilevel++){
        mpTAssyModel= mpGMGModel->getAssyModel(ilevel);

        CContactMesh *pContactMesh= new CContactMesh;// 接合メッシュの生成
        pContactMesh->setID(contactID);
        pContactMesh->setLevel(ilevel);

        pContactMesh->setRank(myRank);
        pContactMesh->setTransmitRank(transRank);

        pContactMesh->setProp(nProp);
        
        mpTAssyModel->addContactMesh(pContactMesh, contactID);
    };
}
// コンタクトノードの生成(Level==0)
//
void CMeshFactory::GeneContactNode(const uint& mgLevel, const uint& contactID, const uint& conNodeID, const vdouble& vCoord,
        const string& s_param_type, const uint& numOfVector, const uint& numOfScalar,
        bool bmesh, const uint& meshID, const uint& nodeID,
        const uint& rank, const uint& maslave)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);

    CContactMesh *pConMesh= mpTAssyModel->getContactMesh_ID(contactID);
    
    CContactNode *pConNode= new CContactNode;
    pConNode->setLevel(mgLevel);
    pConNode->pushLevelMarking();//2010.05.27
    pConNode->setID(conNodeID);
    pConNode->setCoord(vCoord);
    if(bmesh){ pConNode->setMeshID(meshID); pConNode->markingSelfMesh();}
    if(bmesh){ pConNode->setNodeID(nodeID); pConNode->markingSelfNode();}
    pConNode->setRank(rank);

    if(s_param_type=="v" || s_param_type=="V" || s_param_type=="sv" || s_param_type=="SV"){
        pConNode->resizeDisp(numOfVector);
        pConNode->initDisp();
    }
    if(s_param_type=="s" || s_param_type=="S" || s_param_type=="sv" || s_param_type=="SV"){
        pConNode->resizeScalar(numOfVector);
        pConNode->initScalar();
    }

    pConMesh->addConNode(pConNode, conNodeID);

    switch(maslave){
    case(0)://マスターConNode
        pConMesh->addMasterConNode(pConNode,conNodeID);
        break;
    case(1)://スレーブConNode
        pConMesh->addSlaveConNode(pConNode,conNodeID);
        break;
    default:
        break;
    }
}

// マスター面の生成(Level==0)
//
//
void CMeshFactory::GeneMasterFace(const uint& contactID, const uint& shapeType, const uint& masterFaceID,
        bool bmesh, const uint& meshID, const uint& elemID, const uint& elemFaceID,
        const vuint& vConNodeID, const uint& face_rank)
{
    // レベル0の,マスター&スレーブ要素をマーキング
    //  -> レベル0以外は,dividHexa()等でマーキング
    //
    uint mgLevel(0);
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    CSkinFace *pMFace= new CMasterFace;

    //初期化
    pMFace->setShapeType(shapeType);

    // bmesh がtrueの場合,自身のMeshの表面がSkinFaceになっているので,マーキングする.
    if(bmesh){
        CMesh *pMMesh;
        pMMesh= mpTAssyModel->getMesh_ID(meshID);
        
        CElement *pMasterElem;
        pMasterElem= pMMesh->getElement(elemID);
        pMasterElem->markingMPCMaster();
        pMasterElem->markingMPCFace(elemFaceID);
        
        pMFace->markingSelf();//自身のMeshの表面メッシュであることをマーキング
        pMFace->setMeshID(meshID);
        pMFace->setElementID(elemID);
        pMFace->setFaceID(elemFaceID);
    }
    
    CContactMesh *pConMesh= mpTAssyModel->getContactMesh_ID(contactID);
    
    pMFace->setID(masterFaceID);
    pMFace->setRank(face_rank);
    pMFace->setLevel(mgLevel);// mgLevel=0 入力レベル

    uint icnode, numOfConNode(vConNodeID.size());
    for(icnode=0; icnode< numOfConNode; icnode++){
        CContactNode *pConNode= pConMesh->getContactNode_ID(vConNodeID[icnode]);
        pMFace->addNode(pConNode);
    };
    //2次要素の場合、辺ノードをセットしておく
    if(pMFace->getOrder()==ElementOrder::Second){
        uint nNumOfEdge= pMFace->getNumOfEdge();
        uint nNumOfVert= pMFace->getNumOfVert();
        for(uint iedge=0; iedge < nNumOfEdge; iedge++){
            CContactNode *pConNode= pMFace->getNode(nNumOfVert + iedge);
            pMFace->setEdgeConNode(pConNode, iedge);
            pMFace->markingEdgeNode(iedge);
        };
    }


    pConMesh->addMasterFace(pMFace);
}
// スレーブ面の生成(Level==0)
//
//
void CMeshFactory::GeneSlaveFace(const uint& contactID, const uint& shapeType, const uint& slaveFaceID,
        bool bmesh, const uint& meshID, const uint& elemID, const uint& elemFaceID,
        const vuint& vConNodeID, const uint& face_rank)
{
    // レベル0の,マスター&スレーブ要素をマーキング
    //  -> レベル0以外は,dividHexa()等でマーキング
    //
    uint mgLevel(0);
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    CSkinFace *pSFace= new CSkinFace;

    //初期化
    pSFace->setShapeType(shapeType);

    // bmesh がtrueの場合,自身のMeshの表面がSkinFaceになっているので,マーキングする.
    if(bmesh){
        CMesh *pSMesh;
        pSMesh= mpTAssyModel->getMesh_ID(meshID);
        
        CElement *pSlaveElem;
        pSlaveElem= pSMesh->getElement(elemID);
        pSlaveElem->markingMPCSlave();
        pSlaveElem->markingMPCFace(elemFaceID);

        pSFace->markingSelf();//自身のMeshの表面メッシュであることをマーキング
        pSFace->setMeshID(meshID);
        pSFace->setElementID(elemID);
        pSFace->setFaceID(elemFaceID);
    }

    CContactMesh *pConMesh= mpTAssyModel->getContactMesh_ID(contactID);
    
    pSFace->setID(slaveFaceID);
    pSFace->setRank(face_rank);
    pSFace->setLevel(mgLevel);// mgLevel=0 入力レベル
    
    uint icnode, numOfConNode(vConNodeID.size());
    for(icnode=0; icnode< numOfConNode; icnode++){
        CContactNode *pConNode= pConMesh->getContactNode_ID(vConNodeID[icnode]);
        pSFace->addNode(pConNode);
    };
    //2次要素の場合、辺ノードをセットしておく
    if(pSFace->getOrder()==ElementOrder::Second){
        uint nNumOfEdge= pSFace->getNumOfEdge();
        uint nNumOfVert= pSFace->getNumOfVert();
        for(uint iedge=0; iedge < nNumOfEdge; iedge++){
            CContactNode *pConNode= pSFace->getNode(nNumOfVert + iedge);
            pSFace->setEdgeConNode(pConNode, iedge);
            pSFace->markingEdgeNode(iedge);
        };
    }

    pConMesh->addSlaveFace(pSFace);
}

// MPC接触Mesh生成(全ての階層Level) : マスター面,スレーブ面のID番号の管理
//  -> Level=0のマスター & スレーブはセット済み.
//  -> Level=0から始めてprogMeshを利用する.
// --
void CMeshFactory::refineContactMesh()
{
    // Meshリファイン後に処理
    // --
    // ContactMeshをリファイン.
    // --
    CAssyModel  *pAssy,*pProgAssy;
    CContactMesh *pConMesh,*pProgConMesh;
    CSkinFace  *pSkinFace;
    vector<CSkinFace*> vProgFace;//refineで生成されるSkinFaceの子供
    uint maslave;                //マスター,スレーブ切り替えINDEX

    uint meshID,elemID;
    CMesh *pMesh;
    CElement *pElem;

    uint faceID;//progFaceのID番号生成用途
    uint maxLayer;//Octreeの最上位レイヤー
    // ----
    // progAssyのContactMeshにRefineしたContactMeshをセットしていくので,ループは < mMGLevel となる.
    // ----
    uint ilevel;
    for(ilevel=0; ilevel< mMGLevel; ilevel++){
        
        pAssy =  mpGMGModel->getAssyModel(ilevel);    //カレントレベル
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);//上段のレベル
        
        
        uint numOfCont= pAssy->getNumOfContactMesh();
        uint icont;
        for(icont=0; icont< numOfCont; icont++){
            pConMesh= pAssy->getContactMesh(icont);
            pProgConMesh= pProgAssy->getContactMesh(icont);
            
            pConMesh->setupCoarseConNode(pProgConMesh);      //現在LevelのConMeshのノードを上位のConMeshに丸ごとセット
            pConMesh->setupAggSkinFace();                    //ContactNode周囲のSkinFaceIDを収集
            pConMesh->setupEdgeConNode(pProgConMesh, ilevel);//辺ノードの生成,辺接続Faceのセット,新ノードをprogConMeshに追加,IDカウント
            pConMesh->setupFaceConNode(pProgConMesh);        //面ノードの生成,新ノードをprogConMeshに追加,IDカウント
            
            uint numOfSkinFace;
            uint iface;
            //マスター,スレーブ切り替えループ
            for(maslave=0; maslave< 2; maslave++){
                faceID=0;//新たな面IDカウンター(mvLevel別,ContactMesh別,マスター&スレーブ別なので,ここで"0"初期化)
                
                if(maslave==0) numOfSkinFace= pConMesh->getNumOfMasterFace();//マスター面数
                if(maslave==1) numOfSkinFace= pConMesh->getNumOfSlaveFace(); //スレーブ面数

                for(iface=0; iface< numOfSkinFace; iface++){
                    
                    if(maslave==0)  pSkinFace= pConMesh->getMasterFace(iface);//マスター面
                    if(maslave==1)  pSkinFace= pConMesh->getSlaveFace(iface); //スレーブ面
                    
                    // 自身のMeshに存在するMPC面であれば, RefinしたProgFaceにelemID,elemFaceID,NodeIDをセットする.
                    //
                    if(pSkinFace->isSelf() && pSkinFace->getNumOfEdge()!=0){
                        meshID= pSkinFace->getMeshID(); elemID= pSkinFace->getElementID();
                        pMesh= pAssy->getMesh_ID(meshID);
                        pElem= pMesh->getElement(elemID);

                        pSkinFace->refine(pElem, faceID);//// <<<<<<<<-- 面のRefineと新FaceIDのカウントアップ
                    }else{
                        pSkinFace->refine(NULL, faceID); //// <<<<<<<<-- 面のRefineと新FaceIDのカウントアップ
                    }
                    
                    vProgFace= pSkinFace->getProgFace();//RefineしたSkinFaceを取得
                    
                    if(maslave==0) pProgConMesh->addMasterFace(vProgFace);//progContactMeshにRefineしたSkinFaceを追加セット
                    if(maslave==1) pProgConMesh->addSlaveFace(vProgFace); //  同上
                    
                };//ifaceループ
            };//maslaveループ(マスター,スレーブ切り替え)
            
        };//icontループ(ContactMesh)
    };//ilevelループ(マルチグリッドLevel)



    // 2次要素対応(最終Levelに辺ノード生成)
    // ----
    pAssy =  mpGMGModel->getAssyModel(mMGLevel);//最終レベル
    uint numOfCont= pAssy->getNumOfContactMesh();
    uint icont;
    for(icont=0; icont< numOfCont; icont++){
        pConMesh= pAssy->getContactMesh(icont);

        pConMesh->setupAggSkinFace();              //ContactNode周囲のSkinFaceIDを収集
        pConMesh->setupEdgeConNode(NULL, mMGLevel);//辺ノードの生成,辺接続Faceのセット,新ノードをprogConMeshに追加,IDカウント
    };



    // ContactMeshに,八分木を生成 
    //
    for(ilevel=0; ilevel < mMGLevel+1; ilevel++){
        pAssy= mpGMGModel->getAssyModel(ilevel);
        uint numOfCont= pAssy->getNumOfContactMesh();
        uint icont;
        for(icont=0; icont < numOfCont; icont++){
            pConMesh= pAssy->getContactMesh(icont);

            uint nRange;
            nRange= pConMesh->getNumOfConNode();
            
            uint nDigitCount(0);
            while(nRange > 100){
                nRange /= 10;
                nDigitCount++;
            };
            //maxLayer= ilevel+1; //八分木レイヤー数(八分木テスト用)
            maxLayer= nDigitCount;//八分木レイヤー数
            
            pConMesh->generateOctree(maxLayer);
        };
    };
}



// CommMesh2 (節点共有型 通信テーブル) のRefine
// ----
void CMeshFactory::refineCommMesh2()
{
    CAssyModel *pAssy,*pProgAssy;
    CMesh *pMesh,*pProgMesh;
    CCommMesh2 *pCommMesh2,*pProgCommMesh2;
    CCommFace *pCommFace;
    
    vector<CCommFace*> mvCommFace;
    CCommFace *pProgCommFace;
    uint countID(0);//progCommFaceのID用(Faceは,新規にIDを割り当てるので,"0"から)
    

    uint ilevel;
    // progMeshが最上位Levelになるまでループ
    // --
    for(ilevel=0; ilevel< mMGLevel; ilevel++){

        pAssy= mpGMGModel->getAssyModel(ilevel);
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);
        
        uint imesh, numOfMesh;
        numOfMesh= pAssy->getNumOfMesh();

        for(imesh=0; imesh< numOfMesh; imesh++){

            pMesh= pAssy->getMesh(imesh);
            pProgMesh= pProgAssy->getMesh(imesh);
            
            uint icomm, numOfComm;
            numOfComm= pMesh->getCommMesh2Size();

            for(icomm=0; icomm< numOfComm; icomm++){
                pCommMesh2= pMesh->getCommMesh2IX(icomm);
                
                pProgCommMesh2 = new CCommMesh2;//// <<<<<<<<< CommMesh2生成

                pProgCommMesh2->setLevel(ilevel+1);
                pProgCommMesh2->setID(pCommMesh2->getID());
                pProgCommMesh2->setRank(pCommMesh2->getRank());
                pProgCommMesh2->setTransmitRank(pCommMesh2->getTrasmitRank());
                
                pProgMesh->setCommMesh2(pProgCommMesh2);////<<<<< 上位のMeshへprogCommMesh2をセット

                //--------------------
                //CommMesh2のRefine準備
                //--------------------
                pCommMesh2->setupCommNode(pProgCommMesh2);//上位へCommNodeをセット
                pCommMesh2->setupAggFace();
                pCommMesh2->setupEdgeCommNode(pProgCommMesh2, ilevel);//上位へ辺のCommNodeをセット
                pCommMesh2->setupFaceCommNode(pProgCommMesh2);//上位へ面のCommNodeをセット
                
                
                //-----------------
                //CommMesh2のRefine
                //-----------------
                uint iface, numOfFace;
                numOfFace= pCommMesh2->getCommFaceSize();
                for(iface=0; iface< numOfFace; iface++){
                    pCommFace= pCommMesh2->getCommFaceIX(iface);
                    
                    //CommFaceが載っている要素
                    uint elemID = pCommFace->getElementID();
                    CElement *pElem= pMesh->getElement(elemID);
                    
                    mvCommFace= pCommFace->refine(pElem);////// <<<<<<<<<<<<<<<<<< リファイン

                    uint ipface,numOfProgFace;
                    numOfProgFace= mvCommFace.size();
                    for(ipface=0; ipface< numOfProgFace; ipface++){
                        pProgCommFace= mvCommFace[ipface];
                        pProgCommFace->setID(countID);

                        pProgCommMesh2->addCommFace(pProgCommFace);///////// 新CommMesh2へ分割したFaceをセット

                        countID++;/////// 新IDカウントアップ: CommFace
                    };
                };//ifaceループ
                
                
                //-----------------------------
                // CommFace:Quad,Triangleの場合
                //-----------------------------
                // 1.面ノード：Faceの要素ID-Entity番号をたどって,
                //    MeshのNodeをセット
                //----
                // 2.辺ノード：Faceの要素IDと
                //     辺両端のノード情報から,MeshのNodeをセット
                //----
                uint elemID,entity_num;
                CElement  *pElem;
                CNode     *pFaceNode;
                CCommNode *pFaceCommNode;
                for(iface=0; iface< numOfFace; iface++){

                    pCommFace= pCommMesh2->getCommFaceIX(iface);
                    
                    // FaceCommNodeへNodeをセット
                    //         Quad,Triangle
                    if(pCommFace->getNumOfEdge() > 2){
                        
                        pFaceCommNode= pCommFace->getFaceCommNode();

                        elemID= pCommFace->getElementID();
                        entity_num= pCommFace->getElementFaceID();

                        pElem= pMesh->getElement(elemID);
                        pFaceNode= pElem->getFaceNode(entity_num);

                        pFaceCommNode->setNode(pFaceNode);

                    }//if(NumOfEdge > 2)

                    // 辺のCommNodeへNodeをセット
                    //
                    PairCommNode pairCommNode;
                    CCommNode *pEdgeCommNode;
                    CNode *pNodeFir, *pNodeSec;
                    CNode *pEdgeNode;
                    uint iedge, numOfEdge;
                    numOfEdge= pCommFace->getNumOfEdge();
                    
                    for(iedge=0; iedge< numOfEdge; iedge++){
                        pairCommNode= pCommFace->getEdgePairCommNode(iedge);

                        pNodeFir= pairCommNode.first->getNode();
                        pNodeSec= pairCommNode.second->getNode();
                        
                        uint edgeIndex;
                        edgeIndex= pElem->getEdgeIndex(pNodeFir,pNodeSec);
                        pEdgeNode= pElem->getEdgeInterNode(edgeIndex);

                        pEdgeCommNode= pCommFace->getEdgeCommNode(iedge);
                        pEdgeCommNode->setNode(pEdgeNode);
                        
                    };//iedgeループ

                    //// 2 次要素の辺ノードをmvCommNodeへ移し替え => CommMesh2:setupEdgeCommNode に移動
                    //// pCommFace->replaceEdgeCommNode();
                    
                };//ifaceループ
            };//icomm ループ
        };//imesh ループ
    };// iLevel ループ

    //
    // 2次要素の場合、辺ノードを最終Levelに追加
    // * 1次要素も辺に生成されるが、mvCommNodeには追加されない
    //
    pAssy= mpGMGModel->getAssyModel(mMGLevel);//最終LevelのAssyModel
    uint imesh;
    uint nNumOfMesh= pAssy->getNumOfMesh();
    for(imesh=0; imesh < nNumOfMesh; imesh++){
        pMesh= pAssy->getMesh(imesh);

        uint icomm;
        uint nNumOfComm= pMesh->getCommMesh2Size();

        for(icomm=0; icomm< nNumOfComm; icomm++){
            pCommMesh2= pMesh->getCommMesh2IX(icomm);

            pCommMesh2->setupAggFace();
            pCommMesh2->setupEdgeCommNode(NULL, ilevel);//上位へ辺のCommNodeをセット

            // ---
            // 辺CommNodeへNodeをセット
            // ---
            // Face
            uint iface;
            uint nNumOfFace = pCommMesh2->getCommFaceSize();
            for(iface=0; iface< nNumOfFace; iface++){
                pCommFace= pCommMesh2->getCommFaceIX(iface);

                //CommFaceが載っている要素
                uint elemID = pCommFace->getElementID();
                CElement *pElem= pMesh->getElement(elemID);

                // 辺のCommNodeへNodeをセット
                //
                PairCommNode pairCommNode;
                CCommNode *pEdgeCommNode;
                CNode *pNodeFir, *pNodeSec;
                CNode *pEdgeNode;
                uint iedge, numOfEdge;
                numOfEdge= pCommFace->getNumOfEdge();

                // Edge
                for(iedge=0; iedge< numOfEdge; iedge++){
                    pairCommNode= pCommFace->getEdgePairCommNode(iedge);

                    pNodeFir= pairCommNode.first->getNode();
                    pNodeSec= pairCommNode.second->getNode();

                    uint edgeIndex;
                    edgeIndex= pElem->getEdgeIndex(pNodeFir,pNodeSec);
                    pEdgeNode= pElem->getEdgeInterNode(edgeIndex);

                    pEdgeCommNode= pCommFace->getEdgeCommNode(iedge);
                    pEdgeCommNode->setNode(pEdgeNode);

                };//iedge loop
            };//iface loop

        };//icomm loop
        
    };//imesh loop
}






// CommMesh2の生成
// --
void CMeshFactory::GeneCommMesh2(const uint& mgLevel, const uint& mesh_id, const uint& comID, 
        const uint& numOfFace, const uint& numOfCommNode,
        const uint& myRank, const uint& nTransmitRank)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);

    mpTCommMesh2= new CCommMesh2;// CommMesh2の生成

    mpTCommMesh2->setLevel(mgLevel);
    mpTCommMesh2->setID(comID);
    mpTCommMesh2->reserveCommFace(numOfFace);
    mpTCommMesh2->reserveCommNode(numOfCommNode);
    mpTCommMesh2->setRank(myRank);
    mpTCommMesh2->setTransmitRank(nTransmitRank);

    mpTMesh->setCommMesh2(mpTCommMesh2);
}

// CommMesh2用途のCommFace生成
// --
void CMeshFactory::GeneCommFace(const uint& mgLevel, const uint& commeshID, const uint& face_id,
            const uint& mesh_id,const uint elem_id, const uint& elem_ent_num, const uint& elem_type, const vuint& vCommNodeID)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh_ID(mesh_id);

    //要素へのマーキング
    CElement *pElement= mpTMesh->getElement(elem_id);
    pElement->markingCommMesh2();
    pElement->markingCommEntity(elem_ent_num);


    mpTCommMesh2= mpTMesh->getCommMesh2(commeshID);

    CCommFace *pCommFace= new CCommFace;// CommFaceの生成

    pCommFace->setID(face_id);
    pCommFace->setElementID(elem_id);
    pCommFace->setElementFaceID(elem_ent_num);
    pCommFace->setMGLevel(mgLevel);
    
    uint nNumOfVert, nNumOfEdge, nOrder;

    switch(elem_type)
    {
        case(ElementType::Quad):
            nNumOfVert= 4; 
            nNumOfEdge= 4;
            nOrder = ElementOrder::First;
            break;

        case(ElementType::Quad2):
            nNumOfVert= 4;
            nNumOfEdge= 4;
            nOrder = ElementOrder::Second;
            break;

        case(ElementType::Triangle):
            nNumOfVert= 3;
            nNumOfEdge= 3;
            nOrder = ElementOrder::First;
            break;

        case(ElementType::Triangle2):
            nNumOfVert= 3;
            nNumOfEdge= 3;
            nOrder = ElementOrder::Second;
            break;

        case(ElementType::Beam):
            nNumOfVert= 2;
            nNumOfEdge= 1;
            nOrder = ElementOrder::First;
            break;

        case(ElementType::Beam2):
            nNumOfVert= 2;
            nNumOfEdge= 1;
            nOrder = ElementOrder::Second;
            break;
            
        default:
            break;
    }
    pCommFace->initialize(nNumOfVert, nNumOfEdge, nOrder);// CommFace初期化


    uint nNumOfNode = vCommNodeID.size();
    uint i, id;
    CCommNode *pCommNode;
    for(i=0; i< nNumOfNode; i++){
        id= vCommNodeID[i];
        pCommNode= mpTCommMesh2->getCommNode(id);

        pCommFace->setCommNode(i, pCommNode);
    };


    bool bSecond(false);
    if(elem_type==ElementType::Quad2)     bSecond=true;
    if(elem_type==ElementType::Triangle2) bSecond=true;
    if(elem_type==ElementType::Beam2)     bSecond=true;

    if(bSecond){
        for(uint iedge=0; iedge< nNumOfEdge; iedge++){
            id= vCommNodeID[nNumOfVert + iedge];
            pCommNode= mpTCommMesh2->getCommNode(id);
            pCommFace->setEdgeCommNode(pCommNode, iedge);
        };
    }
    mpTCommMesh2->addCommFace(pCommFace);
}

// CommMesh2用途のCommNode生成
// --
void CMeshFactory::GeneCommNodeCM2(const uint& mgLevel, const uint& mesh_id, const uint& node_id,const uint& commeshID,
        const uint& comm_node_id, const vdouble& vCoord)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh_ID(mesh_id);
    mpTCommMesh2= mpTMesh->getCommMesh2(commeshID);
    CNode* pNode= mpTMesh->getNode(node_id);

    CCommNode *pCommNode= new CCommNode;// CommNodeの生成

    pCommNode->setID(comm_node_id);
    pCommNode->setNode(pNode);
    pCommNode->setCoord(vCoord);

    mpTCommMesh2->addCommNode(pCommNode);
}


// --
// グループ
// --
//
// ・GroupObjectの生成
// ・GroupID, GroupNameのセット
//
void CMeshFactory::GeneElemGrpOBJ(const uint& mgLevel, const uint& mesh_id, const vuint& vGrpID, vstring& vGrpName)//ElementGroupの生成
{
    CAssyModel *pAssyModel = mpGMGModel->getAssyModel(mgLevel);
    CMesh *pMesh = pAssyModel->getMesh_ID(mesh_id);

    uint i, nNumOfElemGrp = vGrpID.size();
    for(i=0; i < nNumOfElemGrp; i++){

        CElementGroup *pElemGrp = new CElementGroup;

        pElemGrp->setMesh(pMesh);/// Mesh

        uint nGrpID = vGrpID[i];
        pElemGrp->setID(nGrpID);/// ID

        string sGrpName = vGrpName[i];
        pElemGrp->setName(sGrpName);/// Name

        pMesh->addElemGrp(pElemGrp);
    }
}
//
// ・指定GrpIDへパラメーターをセット
//
void CMeshFactory::setElemID_with_ElemGrp(const uint& mgLevel, const uint& mesh_id, const uint& nGrpID, const vuint& vElemID)
{
    CAssyModel *pAssyModel = mpGMGModel->getAssyModel(mgLevel);
    CMesh *pMesh = pAssyModel->getMesh_ID(mesh_id);

    CElementGroup *pElemGrp = pMesh->getElemGrpID(nGrpID);

    uint i, nNumOfElem=vElemID.size();

    for(i=0; i < nNumOfElem; i++){
        pElemGrp->addElementID(vElemID[i]);
    }
}




