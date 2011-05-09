
#include "ElementType.h"

//
//  SkinFace.cpp
//
//
//
//              2009.10.15
//              2009.10.15
//              k.Takeda
#include "Element.h"
#include <vector>
#include "SkinFace.h"
using namespace pmw;


// construct & destruct
//
CSkinFace::CSkinFace()
{
    mbSelfDom=false;//trueならば自身の計算領域内の表面Face
    
    mvNormalVector.resize(3);//面ベクトル
}
CSkinFace::~CSkinFace()
{
    //std::cout << "~CSkinFace" << std::endl;
}

// Shape(形状タイプ) => 辺ノード,隣接Face配列のresize()
//
void CSkinFace::setShapeType(const uint& shapeType)
{
    Utility::CLogger *pLogger;

    switch(shapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            mShapeType= ElementType::Quad;
            mNumOfEdge= NumberOfEdge::Quad();
            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            mShapeType= ElementType::Triangle;
            mNumOfEdge= NumberOfEdge::Triangle();
            break;
        case(ElementType::Beam):case(ElementType::Beam2):
            mShapeType= ElementType::Beam;
            mNumOfEdge= NumberOfEdge::Beam();
        case(ElementType::Point):
            mShapeType= ElementType::Point;
            mNumOfEdge= 0;
        default:
            //Errorメッセージ出力
            pLogger= Utility::CLogger::Instance();
            pLogger->Info(Utility::LoggerMode::Error,"MPC Face, ShapeType Error @CSkinFace::setShapeType");
            break;
    }
    mvEdgeNode.resize(mNumOfEdge);
    mvbEdgeMarking = new bool[mNumOfEdge];
    mvEdgeFace.resize(mNumOfEdge);

    ////debug
    //cout << "construct SkinFace, number of edge = " << mNumOfEdge << endl;

    uint iedge;
    for(iedge=0; iedge< mNumOfEdge; iedge++){
        mvbEdgeMarking[iedge]=false;
    };

}

// タイプ1
// Refineの為の辺ノードと隣接SkinFace 
// --
// 直接,辺番号を指定してセットする
//
void CSkinFace::setEdgeFace(CSkinFace* pFace, const uint& iedge)
{
    mvEdgeFace[iedge]= pFace;
}
void CSkinFace::setEdgeConNode(CContactNode* pEdgeConNode, const uint& iedge)
{
    mvEdgeNode[iedge]= pEdgeConNode;
    mvbEdgeMarking[iedge]= true;
}

// タイプ2
// Refineの為の辺ノードと隣接SkinFace 
// 両端のContactNodeを指定してセットする== 隣接するSkinFaceからの両端ノード情報から,セットする辺を特定
// --
uint& CSkinFace::getEdgeIndex(PairConNode& pairConNode)
{
    //pConNode0,pConNode1のインデックス番号から辺番号を特定
    uint localNum0, localNum1;
    uint numOfNode= mvConNode.size();
    uint icnode;
    for(icnode=0; icnode< numOfNode; icnode++){
        if(mvConNode[icnode]->getID() == pairConNode.first->getID()) localNum0=icnode;
        if(mvConNode[icnode]->getID() == pairConNode.second->getID()) localNum1=icnode;
    };

    
    CEdgeTree* pEdgeTree= CEdgeTree::Instance();// 局所番号の組み合わせで辺番号が特定
    
    switch(mShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            return pEdgeTree->getQuadEdgeIndex(localNum0, localNum1);

        case(ElementType::Triangle):case(ElementType::Triangle2):
            return pEdgeTree->getTriangleEdgeIndex(localNum0, localNum1);

        case(ElementType::Beam):case(ElementType::Beam2):
            return pEdgeTree->getBeamEdgeIndex(localNum0, localNum1);
            
        default:
            return pEdgeTree->getDisagTypeEdgeIndex(localNum0, localNum1);//999が返る.
    }
}

void CSkinFace::setEdgeFace(CSkinFace* pFace, PairConNode& pairConNode)
{
    uint edgeNum= getEdgeIndex(pairConNode);
    
    mvEdgeFace[edgeNum]= pFace;
}
void CSkinFace::setEdgeConNode(CContactNode* pEdgeConNode, PairConNode& pairConNode)
{
    uint edgeNum= getEdgeIndex(pairConNode);
    
    mvEdgeNode[edgeNum]= pEdgeConNode;
}
bool CSkinFace::isEdgeNodeMarking(PairConNode& pairConNode)
{
    uint edgeNum= getEdgeIndex(pairConNode);
    return mvbEdgeMarking[edgeNum];
}
void CSkinFace::markingEdgeNode(PairConNode& pairConNode)
{
    uint edgeNum= getEdgeIndex(pairConNode);

    mvbEdgeMarking[edgeNum]= true;
}


PairConNode CSkinFace::getEdgePairNode(const uint& iedge)
{
    PairConNode pairConNode;
    Utility::CLogger *pLogger;

    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uint* vertNum;
    switch(mShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            vertNum= pEdgeTree->getQuadLocalNodeNum(iedge);
            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            vertNum= pEdgeTree->getTriangleLocalNodeNum(iedge);
            break;
        case(ElementType::Beam):case(ElementType::Beam2):
            vertNum= pEdgeTree->getBeamLocalNodeNum(iedge);
            break;
        default:
            //Errorメッセージ出力
            pLogger= Utility::CLogger::Instance();
            pLogger->Info(Utility::LoggerMode::Error,"MPC Face, ShapeType Error @CSkinFace::getEdgePairNode");
            break;
    }
    
    pairConNode.first = mvConNode[vertNum[0]];
    pairConNode.second= mvConNode[vertNum[1]];
    
    return pairConNode;
}


// ContactNodeを2個指定することによる,辺のConNode 提供
//
CContactNode* CSkinFace::getEdgeConNode(PairConNode& pairConNode)
{
    uint edgeNum= getEdgeIndex(pairConNode);

    return mvEdgeNode[edgeNum];
}




// virtual 自分自身の生成
//  -> refine()でインスタンス生成する際に使用
//
CSkinFace* CSkinFace::generateFace()
{
    mpOtherFace = new CSkinFace;
    return mpOtherFace;
}

// progFaceへ,ElementのNodeIDをセットする.
//  -> refine()で要素の面IDをセットした後で呼び出す.
//
void CSkinFace::setupNodeID_progFace(CElement* pElem, const uint& numOfConNode)
{
    uint ivert,nvert;
    uint vnodeID[2];
    uint elemEdgeIndex;
    CNode *pEdgeNode, *pFaceNode;
    
    ////cout << "SkinFace::setupNodeID_progFace,  aaaaaaaaaaa" << ", ElemFaceID= " << mElementFaceID << endl;

    // *面中心ノードへMeshのNodeIDをセット
    pFaceNode= pElem->getFaceNode(mElementFaceID);//要素面中心のNode*
    
    ////cout << "SkinFace::setupNodeID_progFace,  bbbbbbbbbbbb" << endl;
    
    mpFaceNode->setNodeID(pFaceNode->getID());

    

    // *辺中間ノードへMeshのNodeIDをセット
    // -----------------------------------------------------
    //  自身のConNodeが所有している,頂点のNodeIDを取得
    //   -> Elemの辺番号を取得
    //    -> 辺の中間Nodeを取得
    //     -> 辺中間NodeのIDを取得して自身の辺中間ConNodeにセット
    // -----------------------------------------------------
    for(ivert=0; ivert< numOfConNode; ivert++){
        if(ivert==numOfConNode-1){  nvert=0;}else{  nvert=ivert+1;}
        //自身の頂点二つから,NodeIDを取得
        vnodeID[0]= mvConNode[ivert]->getNodeID(); vnodeID[1]= mvConNode[nvert]->getNodeID();

        elemEdgeIndex= pElem->getEdgeIndex(vnodeID[0],vnodeID[1]);//NodeIDからElemの辺番号を取得
        pEdgeNode= pElem->getEdgeInterNode(elemEdgeIndex);

        //辺のConNodeに,要素の辺NodeのIDをセット{ *辺番号は,始点の頂点番号と同じ }
        mvEdgeNode[ivert]->setNodeID(pEdgeNode->getID());
    };//ivertループ
    
}

// Faceの"再分割" -> SkinFace(スレーブ面)を分割して再生産
// 
void CSkinFace::refine(CElement *pElem, uint& faceID)
{
    CElement *pProgElem;//要素から分割された子要素
    uint nodeID;//要素のNodeID
    uint localNum;//要素の局所ノード番号
    
    uint numOfConNode;//setupNodeID_progFace()用途
    CNode *pEdgeNode; //BeamのNodeIDセットアップ用途
    
    uint iprog;
    CSkinFace* pFace;
    switch(mShapeType){
        //---
        // 四辺形形状
        //---
        case(ElementType::Quad):case(ElementType::Quad2):
            mvProgFace.reserve(4);
            for(iprog=0; iprog< 4; iprog++){
                pFace = generateFace();//new CSkinFace, new MasterFace
                pFace->setRank(mRank);//新Faceは,rankに変化なし:節点界面での領域分割
                pFace->setLevel(mLevel+1);//新FaceはLevelが一段上がる

                pFace->resizeNode(4);
                pFace->setShapeType(ElementType::Quad);//4個の四辺形
                pFace->setID(faceID);//再分割されたFaceのIDセット
                faceID++;//FaceのIDカウントアップ <-- 次のFaceのID
                
                //自身のMeshに存在するSkinFaceの場合の処理(マーキング,MeshID,ElemID,ElemFaceID)
                if(mbSelfDom){
                    pFace->markingSelf();     //計算領域に存在している要素面のSkinFace
                    pFace->setMeshID(mMeshID);//MeshIDは,refineしても変化しない.
                    
                    // ConNodeのNodeIDに一致する頂点にぶら下がる子要素を取得し,ElementIDをセット
                    // --
                    // iprog:子Faceの生成インデックスは,
                    //   * 次の処理で,Faceの頂点番号に合わせてConNodeをセットするので,頂点番号と同義
                    // --
                    nodeID= mvConNode[iprog]->getNodeID();       //iprogは,頂点番号と同義,理由↑
                    pProgElem= pElem->getProgElem_NodeID(nodeID);//MW3-Refineでは,ProgElemは,頂点にぶら下がっている.

                    //ElementIDをセット(要素番号)
                    pFace->setElementID(pProgElem->getID());//progElemのElementIDをセット
                    
                    
                    //ElementFaceIDをセット(面番号)
                    switch(pElem->getType()){
                        case(ElementType::Hexa):case(ElementType::Quad):
                        case(ElementType::Hexa2):case(ElementType::Quad2):
                            pFace->setFaceID(mElementFaceID);//Hexa,Quadは,分割しても面番号は変わらない.
                            break;
                        case(ElementType::Prism):case(ElementType::Prism2):
                            //Prismの四辺形の面の場合(面番号==2,3,4)
                            localNum= pElem->getLocalVertNum(nodeID);//NodeIDに対応する,局所番号を取得
                            switch(localNum){
                                case(0):
                                    if(mElementFaceID==2) pFace->setFaceID(2);
                                    if(mElementFaceID==4) pFace->setFaceID(4);
                                    break;
                                case(1):
                                    if(mElementFaceID==2) pFace->setFaceID(2);
                                    if(mElementFaceID==3) pFace->setFaceID(5);
                                    break;
                                case(2):
                                    if(mElementFaceID==3) pFace->setFaceID(3);
                                    if(mElementFaceID==4) pFace->setFaceID(4);
                                    break;
                                case(3):
                                    if(mElementFaceID==2) pFace->setFaceID(2);
                                    if(mElementFaceID==4) pFace->setFaceID(4);
                                    break;
                                case(4):
                                    if(mElementFaceID==2) pFace->setFaceID(2);
                                    if(mElementFaceID==3) pFace->setFaceID(5);
                                    break;
                                case(5):
                                    if(mElementFaceID==3) pFace->setFaceID(3);
                                    if(mElementFaceID==4) pFace->setFaceID(4);
                                    break;
                            }//switch(頂点番号)
                            break;
                    }//switch(要素の型)
                    
                    ///cout << "SkinFace::refine, aaaaaaaa" << endl;

                    numOfConNode=4;
                    setupNodeID_progFace(pElem, numOfConNode);//辺と面のConNodeへ,メッシュのNodeIDをセットする.
                    
                    ///cout << "SkinFace::refine, bbbbbbbb" << endl;
                }//if(mbSelfDom)
                mvProgFace.push_back(pFace);
            };//iprogループ
            
            
            // 面の頂点番号に合わせて,ProgFaceの配列を決定.
            // ----
            // 上記の要素IDの設定処理と関連しているので変更不可
            // ----
            // Face 0 (vert==0 にぶら下がるSkinFace)
            mvProgFace[0]->setNode(mvConNode[0],0); mvProgFace[0]->setNode(mvEdgeNode[0],1);
            mvProgFace[0]->setNode(mpFaceNode,  2); mvProgFace[0]->setNode(mvEdgeNode[3],3);
            // Face 1 (vert==1 にぶら下がるSkinFace)
            mvProgFace[1]->setNode(mvConNode[1],0); mvProgFace[1]->setNode(mvEdgeNode[1],1);
            mvProgFace[1]->setNode(mpFaceNode,  2); mvProgFace[1]->setNode(mvEdgeNode[0],3);
            // Face 2 (vert==2 にぶら下がるSkinFace)
            mvProgFace[2]->setNode(mvConNode[2],0); mvProgFace[2]->setNode(mvEdgeNode[2],1);
            mvProgFace[2]->setNode(mpFaceNode,  2); mvProgFace[2]->setNode(mvEdgeNode[1],3);
            // Face 3 (vert==3 にぶら下がるSkinFace)
            mvProgFace[3]->setNode(mvConNode[3],0); mvProgFace[3]->setNode(mvEdgeNode[3],1);
            mvProgFace[3]->setNode(mpFaceNode,  2); mvProgFace[3]->setNode(mvEdgeNode[2],3);

            break;
        //---
        // 三角形形状
        //---
        case(ElementType::Triangle):case(ElementType::Triangle2):
            mvProgFace.reserve(3);
            for(iprog=0; iprog< 3; iprog++){

                pFace = generateFace();//new CSkinFace, new MasterFace

                pFace->setRank(mRank);//新Faceは,rankに変化なし:節点界面での領域分割
                pFace->setLevel(mLevel+1);//新Faceは,Level一段あげる

                pFace->resizeNode(4);//四辺形 2010.7.14改修
                pFace->setShapeType(ElementType::Quad);//3個の四辺形
                pFace->setID(faceID);//再分割されたFaceのIDセット
                faceID++;//FaceのIDカウントアップ <-- 次のFaceのID

                if(mbSelfDom){
                    pFace->markingSelf();//計算領域に存在している要素面のSkinFace
                    pFace->setMeshID(mMeshID);//MeshIDは,refineしても変化しない.

                    // ConNodeのNodeIDに一致する頂点の子要素を取得し,ElementIDをセット
                    // --
                    // iprog:子Faceの生成インデックスは,
                    //   * 次の処理で,Faceの頂点番号に合わせてConNodeをセットするので,頂点番号と同義
                    // --
                    nodeID= mvConNode[iprog]->getNodeID(); //iprogは,頂点番号と同義
                    pProgElem= pElem->getProgElem_NodeID(nodeID);//MW3-Refineでは,ProgElemは,頂点にぶら下がっている.

                    //ElementIDをセット
                    pFace->setElementID(pProgElem->getID());//progElemのElementIDをセット

                    //ElementFaceIDをセット(面番号)
                    switch(pElem->getType()){
                        case(ElementType::Tetra):case(ElementType::Tetra2):
                            localNum= pElem->getLocalVertNum(nodeID);//NodeIDに対応する,局所番号を取得
                            switch(localNum){
                                case(0):
                                    if(mElementFaceID==0) pFace->setFaceID(0);
                                    if(mElementFaceID==1) pFace->setFaceID(2);
                                    if(mElementFaceID==3) pFace->setFaceID(4);
                                    break;
                                case(1):
                                    if(mElementFaceID==0) pFace->setFaceID(0);
                                    if(mElementFaceID==1) pFace->setFaceID(2);
                                    if(mElementFaceID==2) pFace->setFaceID(5);
                                    break;
                                case(2):
                                    if(mElementFaceID==0) pFace->setFaceID(0);
                                    if(mElementFaceID==2) pFace->setFaceID(3);
                                    if(mElementFaceID==3) pFace->setFaceID(4);
                                    break;
                                case(3):
                                    if(mElementFaceID==1) pFace->setFaceID(2);
                                    if(mElementFaceID==2) pFace->setFaceID(1);
                                    if(mElementFaceID==3) pFace->setFaceID(4);
                                    break;
                            }//switch(頂点番号)
                            break;
                        case(ElementType::Triangle):case(ElementType::Triangle2):
                            pFace->setFaceID(mElementFaceID);//三角形は,面番号=0 のみなので,面番号に変化なし.
                            break;
                        case(ElementType::Prism):case(ElementType::Prism2):
                            //Prismの三角形の面の場合(面番号==0,1)
                            // * 頂点 0,1,2 => 面番号 0
                            // * 頂点 3,4,5 => 面番号 1
                            localNum= pElem->getLocalVertNum(nodeID);//NodeIDに対応する,局所番号を取得
                            switch(localNum){
                                case(0):case(1):case(2):
                                    pFace->setFaceID(0);
                                    break;
                                case(3):case(4):case(5):
                                    pFace->setFaceID(1);
                                    break;
                            }//switch(頂点番号)
                            break;
                    }//switch(要素の型)
                    
                    numOfConNode=3;
                    setupNodeID_progFace(pElem, numOfConNode);//辺と面のConNodeへ,メッシュのNodeIDをセットする.
                    
                }//if(mbSelfDom)
                mvProgFace.push_back(pFace);
            };//iprogループ
            
            // Face 0 (vert==0 にぶら下がるSkinFace)
            mvProgFace[0]->setNode(mvConNode[0],0); mvProgFace[0]->setNode(mvEdgeNode[0],1);
            mvProgFace[0]->setNode(mpFaceNode,  2); mvProgFace[0]->setNode(mvEdgeNode[2],3);
            // Face 1 (vert==1 にぶら下がるSkinFace)
            mvProgFace[1]->setNode(mvConNode[1],0); mvProgFace[1]->setNode(mvEdgeNode[1],1);
            mvProgFace[1]->setNode(mpFaceNode,  2); mvProgFace[1]->setNode(mvEdgeNode[0],3);
            // Face 2 (vert==2 にぶら下がるSkinFace)
            mvProgFace[2]->setNode(mvConNode[2],0); mvProgFace[2]->setNode(mvEdgeNode[2],1);
            mvProgFace[2]->setNode(mpFaceNode,  2); mvProgFace[2]->setNode(mvEdgeNode[1],3);

            //debug
            //cout << "SkinFace::refine, mvEdgeNode[2] id = " << mvEdgeNode[2]->getID() << endl;

            break;
        //---
        // ビーム形状
        //---
        case(ElementType::Beam):case(ElementType::Beam2):
            mvProgFace.reserve(2);
            for(iprog=0; iprog< 2; iprog++){
                pFace = generateFace();//new CSkinFace, new MasterFace
                pFace->setRank(mRank);//新Faceは,rankに変化なし:節点界面での領域分割
                pFace->setLevel(mLevel+1);
                
                pFace->resizeNode(2);
                pFace->setShapeType(ElementType::Beam);//2個の線分
                pFace->setID(faceID);//再分割されたFaceのIDセット
                faceID++;//FaceのIDカウントアップ <-- 次のFaceのID
                if(mbSelfDom){
                    pFace->markingSelf();//計算領域に存在している要素面のSkinFace
                    pFace->setMeshID(mMeshID);

                    nodeID= mvConNode[iprog]->getNodeID(); //iprogは,頂点番号と同義
                    pProgElem= pElem->getProgElem_NodeID(nodeID);

                    //ElementIDをセット
                    pFace->setElementID(pProgElem->getID());

                    //---
                    //ElementFaceID(面番号)をセット...Beamなので辺番号==0を渡す
                    //---
                    pFace->setFaceID(0);

                    //--
                    //辺のConNodeへ,メッシュのNodeIDをセットする.
                    //--
                    pEdgeNode= pElem->getEdgeInterNode(0);
                    mvEdgeNode[0]->setNodeID(pEdgeNode->getID());
                    mpFaceNode->setNodeID(pEdgeNode->getID());//Beamなので辺中間ノードと同じ

                }//if(mbSelfDom)
                mvProgFace.push_back(pFace);
            };//iprogループ
            
            // Face 0 (Beam 0) (vert==0 にぶら下がるSkinFace)
            mvProgFace[0]->setNode(mvConNode[0],0); mvProgFace[0]->setNode(mvEdgeNode[0],1);
            // Face 1 (Beam 1) (vert==1 にぶら下がるSkinFace)
            mvProgFace[1]->setNode(mvConNode[1],0); mvProgFace[1]->setNode(mvEdgeNode[0],1);

            break;
    }//switch(mShapeType)



    //    //debug
    //    cout << "SkinFace::refine, EdgeConNode :"
    //         <<  mvEdgeNode[0]->getID() << ", "
    //         <<  mvEdgeNode[1]->getID() << ", "
    //         <<  mvEdgeNode[2]->getID() << ", "
    //         <<  mvEdgeNode[3]->getID() << endl;
    //    //debug
    //    cout << "SkinFace::refine, progFace : " << endl;
    //    if(mShapeType==ElementType::Quad){
    //        uint i;
    //        for(i=0; i < 4; i++)
    //            cout << ", face ID= " << mvProgFace[i]->getID()
    //                 << ", conID= " << mvProgFace[i]->getNode(0)->getID()
    //                 << ", conID= " << mvProgFace[i]->getNode(1)->getID()
    //                 << ", conID= " << mvProgFace[i]->getNode(2)->getID()
    //                 << ", conID= " << mvProgFace[i]->getNode(3)->getID()
    //                 << endl;
    //    }
}


// 面ベクトル(正規化)
// --
// 計算
vdouble& CSkinFace::CalcNzNormalVector()
{
    // 外積に使用する2つのベクトル
    double X1= mvConNode[1]->getX() - mvConNode[0]->getX();
    double Y1= mvConNode[1]->getY() - mvConNode[0]->getY();
    double Z1= mvConNode[1]->getZ() - mvConNode[0]->getZ();

    double X2= mvConNode[2]->getX() - mvConNode[0]->getX();
    double Y2= mvConNode[2]->getY() - mvConNode[0]->getY();
    double Z2= mvConNode[2]->getZ() - mvConNode[0]->getZ();

    // 外積
    // x = y1 z2 - z1 y2
    // y = z1 x2 - x1 z2
    // z = x1 y2 - y1 x2
    //
    mvNormalVector[0]= Y1*Z2 - Z1*Y2;
    mvNormalVector[1]= Z1*X2 - X1*Z2;
    mvNormalVector[2]= X1*Y2 - Y1*X2;

    // 正規化
    double vecLength= sqrt(mvNormalVector[0]*mvNormalVector[0] 
                           + mvNormalVector[1]*mvNormalVector[1]
                             + mvNormalVector[2]*mvNormalVector[2]);

    mvNormalVector[0] /= vecLength;
    mvNormalVector[1] /= vecLength;
    mvNormalVector[2] /= vecLength;

    return mvNormalVector;
}

// 取得のみ
vdouble& CSkinFace::getNzNormalVector()
{
    return mvNormalVector;
}


//Master面用途のvirtual method => CSkinFaceでは無効
//
void CSkinFace::addSlaveNode(CContactNode* pConNode)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();//Loggerのインスタンス
    //警告メッセージ
    pLogger->Info(Utility::LoggerMode::Warn, "invalid method CSkinFace::addSlaveNode");
}
void CSkinFace::CalcSlave(const uint& islave, const uint& valType)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();//Loggerのインスタンス
    //警告メッセージ
    pLogger->Info(Utility::LoggerMode::Warn, "invalid method CSkinFace::CalcSlave");
}
double& CSkinFace::getCoef(const uint& slaveID, const uint& ivert)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();//Loggerのインスタンス
    //警告メッセージ
    pLogger->Info(Utility::LoggerMode::Warn, "invalid method CSkinFace::getCoef");

    return pLogger->getDDummyValue();
}


// 辺ノードvectorのメモリー解放
//
void CSkinFace::deleteProgData()
{
    vector<CContactNode*>().swap(mvEdgeNode);
    vector<CSkinFace*>().swap(mvEdgeFace);

    delete []mvbEdgeMarking;
}













