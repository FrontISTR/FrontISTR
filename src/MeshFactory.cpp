
#include "AssyModel.h"


#include "Element.h"


#include "Mesh.h"


#include "GMGModel.h"

//
//  MeshFactory.cpp
//
//
//
//                      2009.07.23
//			2008.11.05
//			k.Takeda
#include "MeshFactory.h"
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
void CMeshFactory::refineMesh()
{
    // ---
    // 階層数==0であってもCMesh::setupAggregateをコールするので,
    // このルーチン(refineMesh)は,必須のルーチンである.
    // ---

    CAssyModel *pAssy, *pProgAssy;
    CMesh      *pMesh, *pProgMesh;

    CElement         *pElem;//元になった要素
    vector<CElement*> vProgElem;//分割された新しい要素達
    CElement         *pProgElem;//分割された新しい要素

    // 階層”0”のMesh数を基準に各階層のAssyModelのMeshを決める.
    // ---
    pAssy= mpGMGModel->getAssyModel(0);
    uint numOfMesh= pAssy->getNumOfMesh();
    
    uint ilevel,imesh,ielem;
    // ---
    // 階層Level ループ
    // ---
    for(ilevel=0; ilevel< mNumOfMGLevel; ilevel++){
        
        pAssy= mpGMGModel->getAssyModel(ilevel);

        //debug
        cout << "pAssy MultiGrid Level= " << pAssy->getMGLevel() << endl;

        
        // prolongation AssyModel
        //
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);//FileReadRefineブロックでAssyModelは生成済み
        pProgAssy->resizeMesh(numOfMesh);
        // ---
        // Mesh(パーツ) ループ in AssyMode
        // ---
        for(imesh=0; imesh< numOfMesh; imesh++){

            //debug
            cout << "CMeshFactory::refineMesh, ilevel = " << ilevel << endl;
            cout << "CMeshFactory::refineMesh, imesh  = " << imesh  << endl;

            pMesh= pAssy->getMesh(imesh);   //Current_Level Mesh(最初期はLevel==0：ファイル読み込み時)

            //debug
            cout << "pMesh MultiGrid Level= " << pMesh->getMGLevel() << endl;
            
            // CMeshの生成(prolongation Mesh)
            //
            pProgMesh = new CMesh;          //Upper_Level Mesh ==(prolongation Mesh)
            pProgMesh->setMGLevel(ilevel+1);//上位MeshのMultiGridレベルを設定(初期pMeshはファイル読み込み時のLevel==0)
            pProgMesh->setMeshID(pMesh->getMeshID());//Mesh_ID は,同一のIDとする.

            
            // Refineの準備, upLevel節点(上のLevelの節点)生成
            //
            pMesh->setupAggregate(); //Node集合Element, Node集合Nodeの計算
                mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupAggElement finish");
            pMesh->presetProgMesh(pProgMesh);//prolongation_Meshのノード,要素リザーブ(reserve) && pMeshのノード,要素をセット
                mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->presetProgMesh finish");
            pMesh->setupEdgeElement(pProgMesh);//辺(Edge)節点, 辺に集合する要素の計算
                mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupEdgeElement finish");
            pMesh->setupFaceElement(pProgMesh);//面(Face)節点, 面に隣接する要素の計算
                mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupFaceElement finish");
            pMesh->setupVolumeNode(pProgMesh); //体(Volume)節点:要素中心の節点
                mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupVolumeNode finish");
            
            
            uint numOfElem = pMesh->getNumOfElement();
            uint indexCount= 0;// 新たに生成される要素のIDを割振る. -> 各divid関数でカウントアップ.
                               // ElementのID(Index)初期値は,土台のMeshとは無関係 <= Nodeとは異なる.

            for(ielem=0; ielem< numOfElem; ielem++){
                
                pElem= pMesh->getElement(ielem);
                
                vProgElem.clear();// 分割 Elementコンテナのクリア
                
                uint i;
                // divid Element(要素の分割)
                switch(pElem->getType()){
                    case(ElementType::Hexa):
                        
                        vProgElem.reserve(8);//分割された新しい要素の生成
                        for(i=0; i< 8; i++){
                            pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                            vProgElem.push_back(pProgElem);
                        };
                        dividHexa(pElem,vProgElem, indexCount, pProgMesh);
                        
                        break;
                    case(ElementType::Tetra):

                        vProgElem.reserve(4);//分割された新しい要素の生成
                        for(i=0; i< 4; i++){
                            pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                            vProgElem.push_back(pProgElem);
                        };
                        dividTetra(pElem,vProgElem, indexCount, pProgMesh);

                        break;
                    case(ElementType::Prism):

                        vProgElem.reserve(6);//分割された新しい要素の生成
                        for(i=0; i< 6; i++){
                            pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                            vProgElem.push_back(pProgElem);
                        };
                        dividPrism(pElem,vProgElem, indexCount,pProgMesh);

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
                        dividPyramid(pElem,vProgElem, indexCount,pProgMesh);

                        break;
                    case(ElementType::Quad):

                        vProgElem.reserve(4);//分割された新しい要素の生成
                        for(i=0; i< 4; i++){
                            pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                            vProgElem.push_back(pProgElem);
                        };
                        dividQuad(pElem,vProgElem, indexCount,pProgMesh);
                        
                        break;
                    case(ElementType::Triangle):

                        vProgElem.reserve(3);//分割された新しい要素の生成
                        for(i=0; i< 3; i++){
                            pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                            vProgElem.push_back(pProgElem);
                        };
                        dividTriangle(pElem,vProgElem, indexCount,pProgMesh);

                        break;
                    case(ElementType::Beam):

                        vProgElem.reserve(2);//分割された新しい要素の生成
                        for(i=0; i< 2; i++){
                            pProgElem= new CBeam; pProgElem->setMGLevel(ilevel+1);
                            vProgElem.push_back(pProgElem);
                        };
                        dividBeam(pElem,vProgElem, indexCount,pProgMesh);
                        
                        break;
                }//switch エンド

                uint nBaseType;
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
            
            // pProgMeshの
            //   AggregateElementの生成
            //   AggregateNodeの生成
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
            //☆ CMesh::setupAggregate()は,このルーチンの先頭のRefine準備でpMeshに対してコールするのでpProgMeshにはコールしない.
            
            // prolongation AssyModelに,pProgMeshをセット
            //
            pProgAssy->setMesh(pProgMesh,pMesh->getMeshID());//progAssyModel に,progMeshをセット
            
            //debug
            cout << "pProgMesh ノード数 == " << pProgMesh->getNumOfNode() << endl;
            
        };//imesh ループ エンド
    };//ilevel ループ エンド


    // 最終LevelのMeshにNode集合Node,Node集合Elementをセットする.
    //
    pAssy= mpGMGModel->getAssyModel(mNumOfMGLevel);//最終LevelのAssyModel
    
    numOfMesh= pAssy->getNumOfMesh();
    for(imesh=0; imesh< numOfMesh; imesh++){
        pMesh= pAssy->getMesh(imesh);
        pMesh->setupAggregate();
    };
}


// 要素を分割して生成, for refineMesh()
// --
// Hexa(6面体)の分割
//
void CMeshFactory::dividHexa(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;//頂点のノード
    vector<CNode*> vEdgeNode;//辺のノード
    vector<CNode*> vFaceNode;//面のノード
    CNode          *pVolNode;//体中心のノード
    
    //    uint numOfVert,numOfEdge,numOfFace;//各要素の属性(分割の際に使用)
    //    numOfVert= NumberOfVertex::Hexa(); numOfFace= NumberOfFace::Hexa(); numOfEdge= NumberOfEdge::Hexa();

    ////debug
    //cout << "dividHexa start" << endl;

    uint i;
    //頂点のノード
    vVertNode.resize(8); for(i=0; i< 8; i++){ vVertNode[i] = pElem->getNode(i);}
    //辺のノード
    vEdgeNode.resize(12);for(i=0; i< 12; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    //面のノード
    vFaceNode.resize(6); for(i=0; i< 6; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    //体ノード
    pVolNode = pElem->getVolumeNode();

    ////debug
    //cout << "dividHexa inter" << endl;


    // 8個のHexaを生成
    // 要素 0
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    vProgElem[0]->setNode(vFaceNode[0],2); vProgElem[0]->setNode(vEdgeNode[3],3);
    vProgElem[0]->setNode(vEdgeNode[8],4); vProgElem[0]->setNode(vFaceNode[4],5);
    vProgElem[0]->setNode(pVolNode,    6); vProgElem[0]->setNode(vFaceNode[3],7);

    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    vProgElem[1]->setNode(vEdgeNode[1],2); vProgElem[1]->setNode(vFaceNode[0],3);
    vProgElem[1]->setNode(vFaceNode[4],4); vProgElem[1]->setNode(vEdgeNode[9],5);
    vProgElem[1]->setNode(vFaceNode[2],6); vProgElem[1]->setNode(pVolNode,    7);

    // 要素 2
    vProgElem[2]->setNode(vEdgeNode[8],0); vProgElem[2]->setNode(vFaceNode[4],1);
    vProgElem[2]->setNode(pVolNode,    2); vProgElem[2]->setNode(vFaceNode[3],3);
    vProgElem[2]->setNode(vVertNode[4],4); vProgElem[2]->setNode(vEdgeNode[4],5);
    vProgElem[2]->setNode(vFaceNode[1],6); vProgElem[2]->setNode(vEdgeNode[7],7);

    // 要素 3
    vProgElem[3]->setNode(vFaceNode[4],0); vProgElem[3]->setNode(vEdgeNode[9],1);
    vProgElem[3]->setNode(vFaceNode[2],2); vProgElem[3]->setNode(pVolNode,    3);
    vProgElem[3]->setNode(vEdgeNode[4],4); vProgElem[3]->setNode(vVertNode[5],5);
    vProgElem[3]->setNode(vEdgeNode[5],6); vProgElem[3]->setNode(vFaceNode[1],7);

    // 要素 4
    vProgElem[4]->setNode(vEdgeNode[3],0); vProgElem[4]->setNode(vFaceNode[0],1);
    vProgElem[4]->setNode(vEdgeNode[2],2); vProgElem[4]->setNode(vVertNode[3],3);
    vProgElem[4]->setNode(vFaceNode[3],4); vProgElem[4]->setNode(pVolNode,    5);
    vProgElem[4]->setNode(vFaceNode[5],6); vProgElem[4]->setNode(vEdgeNode[11],7);

    // 要素 5
    vProgElem[5]->setNode(vFaceNode[0],0); vProgElem[5]->setNode(vEdgeNode[1],1);
    vProgElem[5]->setNode(vVertNode[2],2); vProgElem[5]->setNode(vEdgeNode[2],3);
    vProgElem[5]->setNode(pVolNode,    4); vProgElem[5]->setNode(vFaceNode[2],5);
    vProgElem[5]->setNode(vEdgeNode[10],6);vProgElem[5]->setNode(vFaceNode[5],7);

    // 要素 6
    vProgElem[6]->setNode(vFaceNode[3],0); vProgElem[6]->setNode(pVolNode,    1);
    vProgElem[6]->setNode(vFaceNode[5],2); vProgElem[6]->setNode(vEdgeNode[11],3);
    vProgElem[6]->setNode(vEdgeNode[7],4); vProgElem[6]->setNode(vFaceNode[1],5);
    vProgElem[6]->setNode(vEdgeNode[6],6); vProgElem[6]->setNode(vVertNode[7],7);

    // 要素 7
    vProgElem[7]->setNode(pVolNode,    0); vProgElem[7]->setNode(vFaceNode[2],1);
    vProgElem[7]->setNode(vEdgeNode[10],2); vProgElem[7]->setNode(vFaceNode[5],3);
    vProgElem[7]->setNode(vFaceNode[1],4); vProgElem[7]->setNode(vEdgeNode[5],5);
    vProgElem[7]->setNode(vVertNode[6],6); vProgElem[7]->setNode(vEdgeNode[6],7);


    // IDのセット
    for(i=0; i< 8; i++){
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット
        vProgElem[i]->setID(indexCount);          // <= indexCountは,数を数えているので,配列Indexは直前の配列数
        ++indexCount;//直線の番号を渡した後で個数を数える

        pProgMesh->setElement(vProgElem[i]);
    };

    ////debug
    //cout << "dividHexa end" << endl;

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
    // 要素 1
    vProgElem[1]->setNode(vVertNode[2],0); vProgElem[1]->setNode(vEdgeNode[2],1);
    vProgElem[1]->setNode(vFaceNode[0],2); vProgElem[1]->setNode(vEdgeNode[1],3);
    vProgElem[1]->setNode(vEdgeNode[5],4); vProgElem[1]->setNode(vFaceNode[3],5);
    vProgElem[1]->setNode(pVolNode,    6); vProgElem[1]->setNode(vFaceNode[2],7);
    // 要素 2
    vProgElem[2]->setNode(vFaceNode[0],0); vProgElem[2]->setNode(vEdgeNode[0],1);
    vProgElem[2]->setNode(vVertNode[1],2); vProgElem[2]->setNode(vEdgeNode[1],3);
    vProgElem[2]->setNode(pVolNode,    4); vProgElem[2]->setNode(vFaceNode[1],5);
    vProgElem[2]->setNode(vEdgeNode[4],6); vProgElem[2]->setNode(vFaceNode[2],7);
    // 要素 3
    vProgElem[3]->setNode(vFaceNode[3],0); vProgElem[3]->setNode(vEdgeNode[3],1);
    vProgElem[3]->setNode(vFaceNode[1],2); vProgElem[3]->setNode(pVolNode,    3);
    vProgElem[3]->setNode(vEdgeNode[5],4); vProgElem[3]->setNode(vVertNode[3],5);
    vProgElem[3]->setNode(vEdgeNode[4],6); vProgElem[3]->setNode(vFaceNode[2],7);

    ////debug
    //cout << "要素にノードをセット@dividTetra" << endl;

    // IDのセット
    for(i=0; i< 4; i++){
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };
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
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[1],0); vProgElem[1]->setNode(vVertNode[0],1);
    vProgElem[1]->setNode(vEdgeNode[0],2); vProgElem[1]->setNode(vFaceNode[0],3);
    vProgElem[1]->setNode(vFaceNode[4],4); vProgElem[1]->setNode(vEdgeNode[3],5);
    vProgElem[1]->setNode(vFaceNode[2],6); vProgElem[1]->setNode(pVolNode,    7);
    // 要素 2
    vProgElem[2]->setNode(vFaceNode[0],0); vProgElem[2]->setNode(vEdgeNode[0],1);
    vProgElem[2]->setNode(vVertNode[1],2); vProgElem[2]->setNode(vEdgeNode[2],3);
    vProgElem[2]->setNode(pVolNode,    4); vProgElem[2]->setNode(vFaceNode[2],5);
    vProgElem[2]->setNode(vEdgeNode[4],6); vProgElem[2]->setNode(vFaceNode[3],7);

    // 要素 3
    vProgElem[3]->setNode(vEdgeNode[5],0); vProgElem[3]->setNode(vFaceNode[4],1);
    vProgElem[3]->setNode(pVolNode,    2); vProgElem[3]->setNode(vFaceNode[3],3);
    vProgElem[3]->setNode(vVertNode[5],4); vProgElem[3]->setNode(vEdgeNode[8],5);
    vProgElem[3]->setNode(vFaceNode[1],6); vProgElem[3]->setNode(vEdgeNode[7],7);
    // 要素 4
    vProgElem[4]->setNode(vFaceNode[4],0); vProgElem[4]->setNode(vEdgeNode[3],1);
    vProgElem[4]->setNode(vFaceNode[2],2); vProgElem[4]->setNode(pVolNode,    3);
    vProgElem[4]->setNode(vEdgeNode[8],4); vProgElem[4]->setNode(vVertNode[3],5);
    vProgElem[4]->setNode(vEdgeNode[6],6); vProgElem[4]->setNode(vFaceNode[1],7);
    // 要素 5
    vProgElem[5]->setNode(pVolNode,    0); vProgElem[5]->setNode(vFaceNode[2],1);
    vProgElem[5]->setNode(vEdgeNode[4],2); vProgElem[5]->setNode(vFaceNode[3],3);
    vProgElem[5]->setNode(vFaceNode[1],4); vProgElem[5]->setNode(vEdgeNode[6],5);
    vProgElem[5]->setNode(vVertNode[4],6); vProgElem[5]->setNode(vEdgeNode[7],7);

    // IDのセット
    for(i=0; i< 6; i++){
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };
}
// ピラミッドの分割
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
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    vProgElem[1]->setNode(vEdgeNode[1],2); vProgElem[1]->setNode(vFaceNode[0],3);
    vProgElem[1]->setNode(vFaceNode[4],4); vProgElem[1]->setNode(vEdgeNode[4],5);
    vProgElem[1]->setNode(vFaceNode[1],6); vProgElem[1]->setNode(pVolNode,    7);
    // 要素 2
    vProgElem[2]->setNode(vFaceNode[0],0); vProgElem[2]->setNode(vEdgeNode[1],1);
    vProgElem[2]->setNode(vVertNode[2],2); vProgElem[2]->setNode(vEdgeNode[2],3);
    vProgElem[2]->setNode(pVolNode,    4); vProgElem[2]->setNode(vFaceNode[1],5);
    vProgElem[2]->setNode(vEdgeNode[5],6); vProgElem[2]->setNode(vFaceNode[2],7);
    // 要素 3
    vProgElem[3]->setNode(vEdgeNode[3],0); vProgElem[3]->setNode(vFaceNode[0],1);
    vProgElem[3]->setNode(vEdgeNode[2],2); vProgElem[3]->setNode(vVertNode[3],3);
    vProgElem[3]->setNode(vFaceNode[3],4); vProgElem[3]->setNode(pVolNode,    5);
    vProgElem[3]->setNode(vFaceNode[2],6); vProgElem[3]->setNode(vEdgeNode[6],7);

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
    // 要素 5 (pyramid)
    vProgElem[5]->setNode(vEdgeNode[4],0); vProgElem[5]->setNode(vVertNode[4],1);
    vProgElem[5]->setNode(vEdgeNode[5],2); vProgElem[5]->setNode(vFaceNode[1],3);
    vProgElem[5]->setNode(pVolNode,    4);
    // 要素 6 (pyramid)
    vProgElem[6]->setNode(vEdgeNode[5],0); vProgElem[6]->setNode(vVertNode[4],1);
    vProgElem[6]->setNode(vEdgeNode[6],2); vProgElem[6]->setNode(vFaceNode[2],3);
    vProgElem[6]->setNode(pVolNode,    4);
    // 要素 7 (pyramid)
    vProgElem[7]->setNode(vEdgeNode[6],0); vProgElem[7]->setNode(vVertNode[4],1);
    vProgElem[7]->setNode(vEdgeNode[7],2); vProgElem[7]->setNode(vFaceNode[3],3);
    vProgElem[7]->setNode(pVolNode,    4);

    // IDのセット
    for(i=0; i< 8; i++){
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット
        
        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };
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
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    vProgElem[1]->setNode(vEdgeNode[1],2); vProgElem[1]->setNode(vFaceNode[0],3);
    // 要素 2
    vProgElem[2]->setNode(vEdgeNode[1],0); vProgElem[2]->setNode(vVertNode[2],1);
    vProgElem[2]->setNode(vEdgeNode[2],2); vProgElem[2]->setNode(vFaceNode[0],3);
    // 要素 3
    vProgElem[3]->setNode(vEdgeNode[2],0); vProgElem[3]->setNode(vVertNode[3],1);
    vProgElem[3]->setNode(vEdgeNode[3],2); vProgElem[3]->setNode(vFaceNode[0],3);

    // IDのセット
    for(i=0; i< 4; i++){
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };
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
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[1],0); vProgElem[1]->setNode(vVertNode[2],1);
    vProgElem[1]->setNode(vEdgeNode[2],2); vProgElem[1]->setNode(vFaceNode[0],3);
    // 要素 2
    vProgElem[2]->setNode(vEdgeNode[2],0); vProgElem[2]->setNode(vVertNode[0],1);
    vProgElem[2]->setNode(vEdgeNode[0],2); vProgElem[2]->setNode(vFaceNode[0],3);

    // IDのセット
    for(i=0; i< 3; i++){
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };
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
    // 要素 1
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);

    // IDのセット
    for(i=0; i< 2; i++){
        vProgElem[i]->setParentID(pElem->getID());//親の要素IDをParentIDにセット

        vProgElem[i]->setID(indexCount);//配列Indexは直前の配列数
        ++indexCount;

        pProgMesh->setElement(vProgElem[i]);
    };
}



// setup to BucketMesh in AssyModel
//
void CMeshFactory::setupBucketMesh(const uint& mgLevel, const uint& num_of_mesh, const uint& maxID, const uint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);

    mpTAssyModel->intializeBucket(maxID, minID);
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
void CMeshFactory::GeneMaterial(const uint& material_id, string& name, vuint& vType, vdouble& vValue)
{
    CMaterial *pMaterial = new CMaterial;

    pMaterial->setID(material_id);
    pMaterial->setName(name);

    for(uint i=0; i< vType.size(); i++) pMaterial->setValue(vType[i],vValue[i]);

    
    mpGMGModel->setMaterial(pMaterial);
}









