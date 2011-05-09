
#include <vector>

#include "FaceTree.h"
#include "EdgeTree.h"
#include "ContactNode.h"
#include "ContactMesh.h"
#include "AssyModel.h"


#include "Element.h"


#include "Mesh.h"

//
//  MeshFactory.cpp
//
//
//
//                      2009.09.29
//			2008.11.05
//			k.Takeda
#include "MeshFactory.h"
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



// Refine for MultiGrid
// --
// current Mesh      => pMesh
// prolongation Mesh => pProgMesh
// 
// current CommMesh  => pCommMesh (通信領域メッシュ)
// prolongation CommMesh => pProgCommMesh
// 
void CMeshFactory::refineMesh()
{
    // --- !注意 ---
    // 階層数==0であってもCMesh::setupAggregateをコールするので,
    // このルーチン(refineMesh)は,必須のルーチンである.
    // -------------

    //-- AseeyModel & Mesh --
    CAssyModel *pAssy, *pProgAssy;
    CMesh      *pMesh, *pProgMesh;
    //-- Element
    CElement         *pElem;//元になった要素
    vector<CElement*> vProgElem;//分割された新しい要素達
    CElement         *pProgElem;//分割された新しい要素

    //-- 通信Mesh
    CCommMesh  *pCommMesh, *pProgCommMesh;
    //-- 通信要素(CommElem)
    CCommElement      *pCommElem;//元になったCommElem(親のCommElem)
    vector<CCommElement*> vProgCommElem;//生成されるprogCommElemのコンテナ
    //--
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
                pMesh->setupAggregate(); //Node集合Element, Node集合Nodeの計算
                mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupAggElement finish at ilevel==0");
            }
            pMesh->presetProgMesh(pProgMesh);//prolongation_Meshのノード,要素リザーブ(reserve) && pMeshのノード,要素をセット
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->presetProgMesh finish");
            pMesh->setupEdgeElement(pProgMesh);//辺(Edge)節点, 辺に集合する要素の計算
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupEdgeElement finish");
            pMesh->setupFaceElement(pProgMesh);//面(Face)節点, 面に隣接する要素の計算
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupFaceElement finish");
            pMesh->setupVolumeNode(pProgMesh); //体(Volume)節点:要素中心の節点
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupVolumeNode finish");
            
            // 新ElementID
            uint numOfElem = pMesh->getNumOfElement();
            uint elementID= 0;// 新たに生成される要素のIDを割振る. -> 各divid関数でカウントアップ.
                               // ElementのID(Index)初期値は,土台のMeshとは無関係 <= Nodeとは異なる.

            for(ielem=0; ielem< numOfElem; ielem++){
                
                pElem= pMesh->getElement(ielem);
                
                vProgElem.clear();// 分割 Elementコンテナのクリア
                //再分割要素の生成
                GeneProgElem(ilevel, pElem, pProgElem, vProgElem, elementID, pProgMesh);

                uint i, nBaseType;
                //生成されたElementの初期化
                for(i=0; i< vProgElem.size(); i++){
                    pProgElem= vProgElem[i];

                    // 局所ノードによる面構成をセットアップ
                    nBaseType= pProgElem->getEntityType();
                    if(nBaseType==BaseElementType::Solid || nBaseType==BaseElementType::Shell)
                        pProgElem->setupFaceCnvNodes();
                };

            };//ielem ループ エンド
            
            // ノード, 数要素数のセット
            pProgMesh->setupNumOfNode();
            pProgMesh->setupNumOfElement();
            
            
            // pProgMeshの AggregateElement, AggregateNodeの生成
            // --
            uint numOfNode= pProgMesh->getNodeSize();
            CAggregateElement *pAggElem;
            CAggregateNode    *pAggNode;
            pProgMesh->reserveAggregate(numOfNode);
            for(uint iAgg=0; iAgg< numOfNode; iAgg++){
                pAggElem = new CAggregateElement;
                pAggNode = new CAggregateNode;

                pProgMesh->setAggElement(pAggElem);
                pProgMesh->setAggNode(pAggNode);
            };
            // !注意 CMesh::setupAggregate()は,このルーチンの先頭のRefine準備でpMeshに対してコールするのでpProgMeshにはコールしない.
            
            // prolongation AssyModelに,pProgMeshをセット
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


            // progCommMeshの前処理
            // --
            pProgMesh->setupAggregate();


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

                ////debug
                //cout << "refineMesh CommMeshループエンド" << endl;

            };//CommMeshループ・エンド
            
            // Mesh のNode,Elementの計算領域整理:MeshのmvNode,mvElementから計算に使用しないNode(DNode),Element(DElement)を移動
            // --
            pProgMesh->sortMesh();

            // Meshが,ソートされたので,Bucketを再セットアップ
            pProgMesh->setupBucketNode();//Bucketの"ID,Index"の一括処理
            pProgMesh->setupBucketElement();

            // <<<< end ::pProgCommCommMeshの生成処理 >>>>
            
        };//imesh ループ エンド
    };//ilevel ループ エンド


    // 最終LevelのMeshにNode集合Node,Node集合Elementをセットする.
    //
    pAssy= mpGMGModel->getAssyModel(mMGLevel);//最終LevelのAssyModel
    
    numOfMesh= pAssy->getNumOfMesh();
    for(imesh=0; imesh< numOfMesh; imesh++){
        pMesh= pAssy->getMesh(imesh);
        pMesh->setupAggregate();
    };
}

// 再分割要素(progElem)の生成
//
void CMeshFactory::GeneProgElem(const uint& ilevel,CElement* pElem, CElement* pProgElem, vector<CElement*>& vProgElem, uint& elementID, CMesh* pProgMesh)
{
    uint i;
    // divid Element(要素の分割)
    switch(pElem->getType()){
        case(ElementType::Hexa):

            vProgElem.reserve(8);//分割された新しい要素の生成
            for(i=0; i< 8; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividHexa(pElem,vProgElem, elementID, pProgMesh);

            break;
        case(ElementType::Tetra):

            vProgElem.reserve(4);//分割された新しい要素の生成
            for(i=0; i< 4; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividTetra(pElem,vProgElem, elementID, pProgMesh);

            break;
        case(ElementType::Prism):

            vProgElem.reserve(6);//分割された新しい要素の生成
            for(i=0; i< 6; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividPrism(pElem,vProgElem, elementID,pProgMesh);

            break;
        case(ElementType::Pyramid):

            vProgElem.reserve(8);//分割された新しい要素の生成
            for(i=0; i< 4; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            for(i=0; i< 4; i++){
                pProgElem= new CPyramid; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividPyramid(pElem,vProgElem, elementID,pProgMesh);

            break;
        case(ElementType::Quad):

            vProgElem.reserve(4);//分割された新しい要素の生成
            for(i=0; i< 4; i++){
                pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividQuad(pElem,vProgElem, elementID,pProgMesh);

            break;
        case(ElementType::Triangle):

            vProgElem.reserve(3);//分割された新しい要素の生成
            for(i=0; i< 3; i++){
                pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividTriangle(pElem,vProgElem, elementID,pProgMesh);

            break;
        case(ElementType::Beam):

            vProgElem.reserve(2);//分割された新しい要素の生成
            for(i=0; i< 2; i++){
                pProgElem= new CBeam; pProgElem->setMGLevel(ilevel+1);
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
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット
        vProgElem[i]->setID(elementID);          // <= elementIDは,いままでの数を数えているので,直前の配列数
        ++elementID;//直前の番号を渡した後でelementIDをカウントアップ

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iface;
    if(pElem->isMPCMaster()){
        for(iface=0; iface< 8; iface++){
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
    if(pElem->isMPCSlave()){
        for(iface=0; iface< 8; iface++){
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
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iface;
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
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iface;
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
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット
        
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
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iprog;
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 4; iprog++){ vProgElem[iprog]->markingMPCMaster(); vProgElem[iprog]->markingMPCFace(0);}
    }
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 4; iprog++){ vProgElem[iprog]->markingMPCSlave(); vProgElem[iprog]->markingMPCFace(0);}
    }
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
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iprog;
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 3; iprog++){ vProgElem[iprog]->markingMPCMaster(); vProgElem[iprog]->markingMPCFace(0);}
    }
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 3; iprog++){ vProgElem[iprog]->markingMPCSlave(); vProgElem[iprog]->markingMPCFace(0);}
    }
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
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };

    // MPC面の属性セット
    uint iprog;
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 2; iprog++) vProgElem[iprog]->markingMPCMaster();
    }
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 2; iprog++) vProgElem[iprog]->markingMPCSlave();
    }
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
void CMeshFactory::GeneMesh(const uint& mgLevel, const uint& mesh_id, const uint& index)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);

    mpTMesh = new CMesh();
    mpTMesh->setMeshID(mesh_id);
    mpTMesh->setMGLevel(mgLevel);

    mpTAssyModel->setBucket(mesh_id, index);

    mpTAssyModel->setMesh(mpTMesh,index);

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
            pNode->reserveScalar(numOfScaParam);
            break;
        case(NodeType::Vector):
            pNode = new CVectorNode();
            pNode->reserveVector(numOfVecParam);
            break;
        case(NodeType::ScalarVector):
            pNode = new CScalarVectorNode();
            pNode->reserveScalar(numOfScaParam);
            pNode->reserveVector(numOfVecParam);
            break;
        default:
            //pNode->InitializeNodeADOF(vParam, num_of_param);
            break;
    }
    pNode->setID(id);// id は、連続したIndex番号にする予定(09.06.23)
    //pNode->setIndex(id);//debug : 本来はIndex番号をセット
    pNode->setCoord(coord);
    pNode->setMGLevel(mgLevel);


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
        case(ElementType::Hexa):
            pElement = new CHexa;
            break;
        case(ElementType::Tetra):
            pElement = new CTetra;
            break;
        case(ElementType::Prism):
            pElement = new CPrism;
            break;
        case(ElementType::Pyramid):
            pElement = new CPyramid;
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
        default:
            mpLogger->Info(Utility::LoggerMode::Error,"Error::GeneElement at Factory");
            break;
    }
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
    pElement->setupFaceCnvNodes();// 局所ノードによる面構成をセットアップ
    
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
void CMeshFactory::reserveAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->reserveAggregate(num_of_node);
}

void CMeshFactory::GeneAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    uint i;
    for(i=0; i< num_of_node; i++){
        CAggregateElement *pAggElem = new CAggregateElement;
        CAggregateNode    *pAggNode = new CAggregateNode;

        mpTMesh->setAggElement(pAggElem);//Node周辺の要素集合
        mpTMesh->setAggNode(pAggNode);
    };
}


//----
// Boundary
//----
// Boundary Node
void CMeshFactory::reserveBoundaryNode(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->reserveBoundaryNode(num_of_bnd);
}

// Boundary Face
void CMeshFactory::reserveBoundaryFace(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    
    mpTMesh->reserveBoundaryFace(num_of_bnd);
}

// Boundary Volume
void CMeshFactory::reserveBoundaryVolume(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    mpTMesh->reserveBoundaryVolume(num_of_bnd);
}

// Boundary Node 生成
void CMeshFactory::GeneBoundaryNode(const uint& mgLevel, const uint& mesh_id, const uint& id,
                                    const uint& dof, const uint& bndType, const vdouble& vVal)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    CBoundaryNode *pBNode = new CBoundaryNode();
    uint i;
    switch(bndType){
        case(BoundaryTypeNode::Accel):
            pBNode->initialize(BoundaryTypeNode::Accel, dof);
            for(i=0; i< dof; i++){ pBNode->setValue(vVal[i],i);}
            break;
        case(BoundaryTypeNode::Velo):
            pBNode->initialize(BoundaryTypeNode::Velo, dof);
            for(i=0; i< dof; i++){ pBNode->setValue(vVal[i],i);}
            break;
        case(BoundaryTypeNode::Load):
            pBNode->initialize(BoundaryTypeNode::Load, dof);
            for(i=0; i< dof; i++){ pBNode->setValue(vVal[i],i);}
            break;
        case(BoundaryTypeNode::Disp):
            pBNode->initialize(BoundaryTypeNode::Disp, dof);
            for(i=0; i< dof; i++){ pBNode->setValue(vVal[i],i);}
            break;
        case(BoundaryTypeNode::Temp):
            pBNode->initialize(BoundaryTypeNode::Temp, dof);
            pBNode->setValue(vVal[0],0);
            break;
        case(BoundaryTypeNode::Thermal_Flux):
            pBNode->initialize(BoundaryTypeNode::Thermal_Flux, dof);
            pBNode->setValue(vVal[0], 0);
            break;
        default:
            break;
    }
    pBNode->setID(id);
    mpTMesh->setBoundaryNode(pBNode);
}

// Boundary Face 生成
void CMeshFactory::GeneBoundaryFace(const uint& mgLevel, const uint& mesh_id, const uint& elem_id, const uint& face_id,
                                    const uint& dof, const uint& bndType, const vdouble& vVal)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    CBoundaryFace *pBFace = new CBoundaryFace();
    uint i;
    switch(bndType){
        case(BoundaryTypeFace::Pressure):
            pBFace->initialize(BoundaryTypeFace::Pressure,dof);
            for(i=0; i< dof; i++){ pBFace->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeFace::Temp):
            pBFace->initialize(BoundaryTypeFace::Temp, dof);
            for(i=0; i< dof; i++){ pBFace->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeFace::Thermal_Flux):
            pBFace->initialize(BoundaryTypeFace::Thermal_Flux, dof);
            for(i=0; i< dof; i++){ pBFace->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeFace::TractionVector):
            pBFace->initialize(BoundaryTypeFace::TractionVector, dof);
            for(i=0; i< dof; i++){ pBFace->setValue(vVal[i], i);}
            break;
        default:
            break;
    }
    pBFace->setElementID(elem_id);
    pBFace->setFaceID(face_id);
    mpTMesh->setBoundaryFace(pBFace);

}

// Boundary Volume 生成
void CMeshFactory::GeneBoundaryVolume(const uint& mgLevel, const uint& mesh_id, const uint& id,
                                     const uint& dof, const uint& bndType, const vdouble& vVal)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);

    CBoundaryVolume* pBVolume = new CBoundaryVolume();
    uint i;
    switch(bndType){
        case(BoundaryTypeVolume::Accel):
            pBVolume->initialize(BoundaryTypeVolume::Accel, dof);
            for(i=0; i< dof; i++){ pBVolume->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeVolume::Centrifugal_Force):
            pBVolume->initialize(BoundaryTypeVolume::Centrifugal_Force, dof);
            for(i=0; i< dof; i++){ pBVolume->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeVolume::Gravity):
            pBVolume->initialize(BoundaryTypeVolume::Gravity, dof);
            for(i=0; i< dof; i++){ pBVolume->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeVolume::Heat):
            pBVolume->initialize(BoundaryTypeVolume::Heat, dof);
            for(i=0; i< dof; i++){ pBVolume->setValue(vVal[i], i);}
            break;
        default:
            break;
    }
    pBVolume->setID(id);
    mpTMesh->setBoundaryVolume(pBVolume);
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
        case(ElementType::Pyramid):
            pCommElem = new CCommPyramid;
            break;
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
        case(ElementType::Pyramid):
            // Pyramid => Hexa縮退 ?
            vProgCommElem.reserve(8);
            for(ivert=0; ivert< 8; ivert++){
                pProgCommElem= new CCommHexa;// <<<<<<<<<-- prolongation CommElement
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);

            break;
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

        ////debug
        //cout << "dividCommElem set ElemID => " << pProgElem->getID() << endl;



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
void CMeshFactory::GeneContactMesh(const uint& contactID)
{
    uint ilevel;
    for(ilevel=0; ilevel< mMGLevel+1; ilevel++){
        mpTAssyModel= mpGMGModel->getAssyModel(ilevel);

        CContactMesh *pContactMesh= new CContactMesh;// 接合メッシュの生成
        pContactMesh->setID(contactID);
        pContactMesh->setLevel(ilevel);
        
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
        const vuint& vConNodeID)
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
    uint icnode, numOfConNode(vConNodeID.size());
    for(icnode=0; icnode< numOfConNode; icnode++){
        CContactNode *pConNode= pConMesh->getContactNode_ID(vConNodeID[icnode]);
        pMFace->addNode(pConNode);
    };

    pConMesh->addMasterFace(pMFace);
}
// スレーブ面の生成(Level==0)
//
//
void CMeshFactory::GeneSlaveFace(const uint& contactID, const uint& shapeType, const uint& slaveFaceID,
        bool bmesh, const uint& meshID, const uint& elemID, const uint& elemFaceID,
        const vuint& vConNodeID)
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
    uint icnode, numOfConNode(vConNodeID.size());
    for(icnode=0; icnode< numOfConNode; icnode++){
        CContactNode *pConNode= pConMesh->getContactNode_ID(vConNodeID[icnode]);
        pSFace->addNode(pConNode);
    };

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
    CAssyModel *pAssy,*pProgAssy;
    CContactMesh *pConMesh,*pProgConMesh;
    uint countID;// <<<<<<<<<<<<<------- 辺,面中心に生成される新たなContactNodeのIDのためのカウンター
    CSkinFace *pSkinFace;
    vector<CSkinFace*> vProgFace;//refineで生成されるSkinFaceの子供
    uint maslave;//マスター,スレーブ切り替えINDEX

    uint meshID,elemID;
    CMesh *pMesh;
    CElement *pElem;

    uint faceID;//progFaceのID番号生成用途
    
    uint ilevel;
    for(ilevel=0; ilevel< mMGLevel; ilevel++){
        pAssy =  mpGMGModel->getAssyModel(ilevel);    //カレントレベル
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);//上段のレベル
        
        //debug
        cout << "mgLevel= " << ilevel << endl;

        uint numOfCont= pAssy->getNumOfContactMesh();
        uint icont;
        for(icont=0; icont< numOfCont; icont++){
            pConMesh= pAssy->getContactMesh(icont);
            pProgConMesh= pProgAssy->getContactMesh(icont);
            
            if(ilevel==0) countID= pConMesh->getNumOfConNode();//新たなノードのID 初期値

            //debug
            cout << "ContactMesh::setupEdgeConNode,setupFaceConNode コール" << endl;

            pConMesh->setupAggSkinFace();//ContactNode周囲のSkinFaceIDを収集
            pConMesh->setupEdgeConNode(pProgConMesh, countID);//辺ノードの生成,辺接続Faceのセット,新ノードをprogConMeshに追加,IDカウント
            pConMesh->setupFaceConNode(pProgConMesh, countID);//面ノードの生成,新ノードをprogConMeshに追加,IDカウント
            
            //debug
            cout << "ContactMesh::setupEdgeConNode,setupFaceConNode   終了" << endl;

            uint numOfSkinFace;
            uint iface;
            
            //マスター,スレーブ切り替えループ
            for(maslave=0; maslave< 2; maslave++){
                faceID=0;//新たな面IDカウンター(mvLevel別,ContactMesh別,マスター&スレーブ別なので,ここで"0"初期化)
                
                if(maslave==0) numOfSkinFace= pConMesh->getNumOfMasterFace();//マスター面数
                if(maslave==1) numOfSkinFace= pConMesh->getNumOfSlaveFace();//スレーブ面数

                for(iface=0; iface< numOfSkinFace; iface++){

                    if(maslave==0)  pSkinFace= pConMesh->getMasterFace(iface);//マスター面
                    if(maslave==1)  pSkinFace= pConMesh->getSlaveFace(iface);//スレーブ面
                    
                    // 自身のMeshに存在するMPC面であれば, RefinしたProgFaceにelemID,elemFaceID,NodeIDをセットする.
                    //
                    if(pSkinFace->isSelf() && pSkinFace->getNumOfEdge()!=0){
                        meshID= pSkinFace->getMeshID(); elemID= pSkinFace->getElementID();
                        pMesh= pAssy->getMesh_ID(meshID);
                        pElem= pMesh->getElement(elemID);

                        pSkinFace->refine(pElem, faceID);//// <<<<<<<<<-- 面のRefineと新FaceIDのカウントアップ
                    }else{
                        pSkinFace->refine(NULL, faceID);//// <<<<<<<<-- 面のRefineと新FaceIDのカウントアップ
                    }
                    vProgFace= pSkinFace->getProgFace();//RefineしたSkinFaceを取得

                    if(maslave==0) pProgConMesh->addMasterFace(vProgFace);//progContactMeshにRefineしたSkinFaceを追加セット
                    if(maslave==1) pProgConMesh->addSlaveFace(vProgFace); //  同上
                    
                };//ifaceループ
            };//maslaveループ(マスター,スレーブ切り替え)
        };//icontループ(ContactMesh)
    };//ilevelループ(マルチグリッドLevel)
}



































