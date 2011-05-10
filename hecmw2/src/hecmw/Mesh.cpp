//
//  Mesh.cpp
//
//                  2009.10.15
//                  2008.11.10
//                  k.Takeda
#include "Vertex.h"
#include "LoggerType.h"
#include "IndexBucket.h"
#include "BoundaryEdgeMesh.h"
#include "BoundaryNodeMesh.h"
#include "BoundaryMesh.h"
#include "CommunicationMesh.h"
#include "Logger.h"
#include "Node.h"

#include <vector>

#include "VectorNode.h"
#include "Element.h"

#include "AggregateElement.h"

#include "Mesh.h"
#include "ScalarNode.h"
#include "ScalarVectorNode.h"
using namespace pmw;

// construct & destruct
// --
CMesh::CMesh(void)
{
    mpLogger = Utility::CLogger::Instance();
    //mpLogger->Info(Utility::LoggerMode::MWDebug, "Use! CMesh(const uint& numOfNode, const uint& numOfElem)");

    //debug
    //mnDummyCount= 0;
}

CMesh::CMesh(const uint& numofnode, const uint& numofelem)
{
    mpLogger = Utility::CLogger::Instance();

    mNumOfNode = numofnode;  mNumOfElement = numofelem;

    mvNode.resize(mNumOfNode);
    mvElement.resize(mNumOfElement);
}

// Node, Element
//
CMesh::~CMesh(void)
{
    //debug
    cout << "~CMesh start,  mMGLevel==" << mMGLevel << endl;

    // 要素集合(AggretateElement)の削除
    //
    if(mnSolutionType==SolutionType::FEM)
        for_each(mvAggElement.begin(), mvAggElement.end(), DeleteObject());

    // ノードの周囲ノード(AggregateNode)の削除
    //
    if(mnSolutionType==SolutionType::FVM)
        for_each(mvAggNode.begin(), mvAggNode.end(), DeleteObject());
    
    // CommMeshの削除
    // --
    for_each(mvCommMesh.begin(), mvCommMesh.end(), DeleteObject());

    // 最大階層Levelのときに全てのNodeを削除
    if(mMGLevel==mMaxMGLevel) for_each(mvNode.begin(), mvNode.end(), DeleteObject());

    // 要素(Element)の削除: '10.08.31
    for_each(mvElement.begin(), mvElement.end(), DeleteObject());


    //debug
    cout << "~CMesh   end,  mMGLevel==" << mMGLevel << endl;
}

// Bucket setup -> Node
//
// case : max_id, min_id 判明
//
void CMesh::initBucketNode(const uint& max_id, const uint& min_id)
{
    moBucket.clearBucketNode();
    moBucket.resizeBucketNode(max_id, min_id);
}
// ID -> Index, setup
void CMesh::setupBucketNodeIndex(const uint& id, const uint& index)
{
    moBucket.setIndexNode(id, index);
}


// case : mvNode setup 済み
//        Max_ID , Min_ID を取得してBucketにセットアップ
//
void CMesh::setupBucketNode()
{
    CNode *pNode;
    uint maxID,minID, i;

    // init
    pNode = mvNode[0];
    maxID = pNode->getID(); minID = pNode->getID();

    // serch Max,Min
    for(i=0; i < mvNode.size(); i++){
        pNode = mvNode[i];
        if(pNode->getID() > maxID) maxID = pNode->getID();
        if(pNode->getID() < minID) minID = pNode->getID();
    };

    moBucket.resizeBucketNode(maxID, minID);

    // set data
    for(i=0; i < mvNode.size(); i++){
        pNode = mvNode[i];
        moBucket.setIndexNode(pNode->getID(), i);
    };
}

// Bucket setup -> Element
//
// case : max_id, min_id 判明
//
void CMesh::initBucketElement(const uint& max_id, const uint& min_id)
{
    moBucket.clearBucketElement();
    moBucket.resizeBucketElement(max_id, min_id);
}
// ID -> Index, setup
void CMesh::setupBucketElementIndex(const uint& id, const uint& index)
{
    moBucket.setIndexElement(id, index);
}


// case : mvElement setup 済み
//        Max_ID , Min_ID を取得してBucketにセットアップ
//
void CMesh::setupBucketElement()
{
    CElement *pElement;
    uint maxID,minID, i;

    // init
    pElement = mvElement[0];
    maxID = pElement->getID(); minID = pElement->getID();

    // serch Max,Min
    for(i=0; i < mvElement.size(); i++){
        pElement = mvElement[i];
        if(pElement->getID() > maxID) maxID = pElement->getID();
        if(pElement->getID() < minID) minID = pElement->getID();
    };

    moBucket.resizeBucketElement(maxID,minID);

    // set data
    for(i=0; i < mvElement.size(); i++){
        pElement = mvElement[i];
        moBucket.setIndexElement(pElement->getID(), i);
    };
}

//--
// Node
//--

// Reserve mvNode
//
void CMesh::reserveNode(const uint& num_of_node)
{
    mvNode.reserve(num_of_node);
    //numOfNode= num_of_node;
}

// Node  setup
// --
// arg::id  ==> Bucket
//
void CMesh::setNode(CNode *pNode)
{
    ////debug
    //cout << "CMesh::setNode  ID = " << pNode->getID() << endl;

    mvNode.push_back(pNode);

    ////debug
    //cout << "CMesh::setNode mvNode.size() = " << mvNode.size() << endl;
}

// IDから配列のIndexを取得して,配列からNode*を取得
//
CNode* CMesh::getNode(const uint &id)
{
    uint index;
    index = moBucket.getIndexNode(id);

    //if(index==-1){ index=id;}//Indexが-1のIDは存在しない.

    return mvNode[index];
}

// 直接Indexを指定してNode*を取得
//
CNode* CMesh::getNodeIX(const uint& index)
{
    return mvNode[index];
}



//--
// Element
//--

// Reserve mvElement
//
void CMesh::reserveElement(const uint& num_of_elem)
{
    mvElement.reserve(num_of_elem);
    //numOfElement= num_of_elem;
}

// Element Setup
// --
// arg::id => Bucket
//
void CMesh::setElement(CElement *pElement)
{
    mvElement.push_back(pElement);
}

// IDからIndex番号を取得して,配列からCElement*を取得
//
CElement* CMesh::getElement(const uint& id)
{
    uint index;
    index = moBucket.getIndexElement(id);

    //if(index==-1){ index = id;}

    return mvElement[index];
}

// 直接Indexを指定してElement*を取得
//
CElement* CMesh::getElementIX(const uint& index)
{
    return mvElement[index];
}



// Aggregate-Element, Aggregate-Node
// ---
//
void CMesh::resizeAggregate(const uint& res_size)
{
//    mvAggElement.resize(res_size + 1);// 配列数とMaxIDがずれることは,ないはずなので＋１は不要…
//    mvAggNode.resize(res_size + 1);   //

    if(mnSolutionType==SolutionType::FEM) mvAggElement.resize(res_size);
    if(mnSolutionType==SolutionType::FVM) mvAggNode.resize(res_size);
}

void CMesh::setAggElement(CAggregateElement* pAggElem, const uint& inode)
{
    mvAggElement[inode] = pAggElem;
}

void CMesh::setAggNode(CAggregateNode* pAggNode, const uint& inode)
{
    mvAggNode[inode] = pAggNode;
}

// Node_ID に対応する,AggregateElementの提供
//
CAggregateElement* CMesh::getAggElem(const uint& node_id)
{
    uint index;
    index = moBucket.getIndexNode(node_id);

    return mvAggElement[index];
}
// 配列Indexに対応する、AggregateElementの提供
//
CAggregateElement* CMesh::getAggElemIX(const uint& inode)
{
    return mvAggElement[inode];
}

// Node_ID に対応する,AggregateNodeの提供
//
CAggregateNode* CMesh::getAggNode(const uint& node_id)
{
    uint index;
    index = moBucket.getIndexNode(node_id);

    return mvAggNode[index];
}
// 配列Indexに対応する、AggregateNodeの提供
//
CAggregateNode* CMesh::getAggNodeIX(const uint& inode)
{
    return mvAggNode[inode];
}

//
// AggregateElement,AggregateNode のセットアップ
//  :nLevelは、2次要素の初期(Level=0)処理のために使用
// 
void CMesh::setupAggregate(const uint& nLevel)
{
    CElement *pElem; CNode *pNode;
    uint ielem, local_id, inode;

    mNumOfNode= mvNode.size();
    mNumOfElement= mvElement.size();


    // 第2段 以降のprolongationのときに,
    //  pNode(Vertex)のElemIndexに同じIDが入らないようにクリア.
    //
    for(inode=0; inode< mNumOfNode; inode++){
        pNode= mvNode[inode];

        pNode->clearAggElemID();
        pNode->clearNeibElemVert();
    };


    // Node(Vertex)が所属する、要素のID番号をセット(in Node) : 2次要素の辺の要素集合については、0段以外は集めない.
    //     (複数の要素に所属している)
    //
    uint numOfLocalNode; uint elemID;
    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem = mvElement[ielem];

        if(nLevel==0 && pElem->getOrder()==ElementOrder::Second){
            // 2次要素( mgLevel == 0 )
            numOfLocalNode= pElem->getNumOfNode();// <== 辺ノードを含めた要素ノード数

            elemID = pElem->getID();

            for(local_id=0; local_id < numOfLocalNode; local_id++){
                pNode = pElem->getNode(local_id);

                pNode->setAggElemID(elemID);
                pNode->setNeibElemVert(elemID, local_id);
            };
        }else{
            // 1次要素 .or. 2次要素( mgLevel >= 1 )
            numOfLocalNode= pElem->getNumOfVert();// <== 頂点数  2010.11.23

            elemID = pElem->getID();

            for(local_id=0; local_id< numOfLocalNode;local_id++){
                pNode = pElem->getNode(local_id);

                pNode->setAggElemID(elemID);      //Node自身に要素ID(AggElementID)をセット
                pNode->setNeibElemVert(elemID,local_id);//Node自身に接続先の頂点番号をセット <= CommMeshのIndex管理で使用
            };
        }
    };

//    // テスト用途
//    // ---------
//    // Node(Vertex)が所属する、要素のID番号をセット(in Node)
//    //
//    //
//    uint numOfLocalNode; uint elemID;
//    for(ielem=0; ielem< mNumOfElement; ielem++){
//        pElem = mvElement[ielem];
//
//        numOfLocalNode= pElem->getNumOfVert();// <== 頂点数  2010.11.23
//
//        elemID = pElem->getID();
//
//        for(local_id=0; local_id< numOfLocalNode;local_id++){
//            pNode = pElem->getNode(local_id);
//
//            pNode->setAggElemID(elemID);      //Node自身に要素ID(AggElementID)をセット
//            pNode->setNeibElemVert(elemID,local_id);//Node自身に接続先の頂点番号をセット <= CommMeshのIndex管理で使用
//        };
//    };


    if(mnSolutionType==SolutionType::FEM){
        CAggregateElement *pAggElem;
        uint iagg_elem, numOfAggElem;
        //
        // Nodeにセットされた要素のIDを基に,AggregateElement(in Mesh)をセット: 2次要素の辺の要素集合については、setupEdgeElement()で集める.
        //
        //for(inode=0; inode< mvNode.size(); inode++){ // △09.08.11にnumOfNodeに変更
        for(inode=0; inode< mNumOfNode; inode++){
            pNode = mvNode[inode];
            numOfAggElem = pNode->getNumOfAggElem();
            
            pAggElem = mvAggElement[inode];// Node Indexに対して,AggElemをセット '10.10.07

            pAggElem->reserve(numOfAggElem);

            for(iagg_elem=0; iagg_elem< numOfAggElem; iagg_elem++){
                elemID = pNode->getAggElemID(iagg_elem);

                pElem= getElement(elemID);

                pAggElem->push(pElem);
            };
        };
    }

    if(mnSolutionType==SolutionType::FVM){
        CAggregateNode *pAggNode;
        vector<CNode*> vConnNode;
        uint ncon, nNumOfConnNode;
        uint iagg_elem, numOfAggElem;
        // VertexのAggregateElementを基に,AggregateNodeをセット
        //
        for(inode=0; inode< mNumOfNode; inode++){
            pNode= mvNode[inode];
            numOfAggElem= pNode->getNumOfAggElem();
            
            pAggNode = mvAggNode[inode];// Node Indexに対して,AggNodeをセット '10.10.07

            pAggNode->reserveNode(numOfAggElem);

            for(iagg_elem=0; iagg_elem< numOfAggElem; iagg_elem++){
                elemID= pNode->getAggElemID(iagg_elem);

                pElem= getElement(elemID);

                // pNodeのIDに対応する辺の一方の端のNodeを拾う.
                vConnNode= pElem->getConnectNode(pNode);
                nNumOfConnNode= vConnNode.size();

                // Node集合にセットしていく.
                //
                for(ncon=0; ncon < nNumOfConnNode; ncon++){
                    pAggNode->setNode(vConnNode[ncon]);//既に取得済みのNodeかどうか,pAggNodeが判定してAggNodeに追加.
                };
            };
        };
    }
}



// ○ 要素辺,面 などにノードを生成する前に,土台Mesh{このMesh自身,つまり"this＊"}の"ノード&要素データ"を上位Meshにプリセット
// ○ pProgMesh::Node,Elementのreserve
//
void CMesh::presetProgMesh(CMesh* pProgMesh)
{
    uint progNumOfNode= mNumOfNode;
    uint progNumOfElem= mNumOfElement;
    //uint numOfHexa(0),numOfTetra(0),numOfPrism(0),numOfPyramid(0),numOfQuad(0),numOfTriangle(0),numOfBeam(0);//現Meshの種類別数
    
    // 増加要素 予測
    CElement *pElem;
    uint ielem;
    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem= mvElement[ielem];
        
        // 要素ごとの要素増加数 & 現Meshの要素種類別数
        switch(pElem->getType()){
            case(ElementType::Hexa):case(ElementType::Hexa2):
                progNumOfElem += 7;//増加する要素数なので、Hexaを8分割できるから、増加数==7;
                //numOfHexa++;
                break;
            case(ElementType::Tetra):case(ElementType::Tetra2):
                progNumOfElem += 3;//増加数==3, 分割数==4
                //numOfTetra++;
                break;
            case(ElementType::Prism):case(ElementType::Prism2):
                progNumOfElem += 5;//増加数==5, 分割数==6
                //numOfPrism++;
                break;
            case(ElementType::Quad):case(ElementType::Quad2):
                progNumOfElem += 3;//増加数==3, 分割数==4
                //numOfQuad++;
                break;
            case(ElementType::Triangle):case(ElementType::Triangle2):
                progNumOfElem += 2;//増加数==2, 分割数==3
                //numOfTriangle++;
                break;
            case(ElementType::Beam):case(ElementType::Beam2):
                progNumOfElem += 1;//増加数==1, 分割数==2
                //numOfBeam++;
                break;
        }
    };

    // 増加ノード予測(大雑把) => 2*2*2=8
    progNumOfNode *= 8;

    pProgMesh->reserveNode(progNumOfNode);   //暫定的なノード数である.
    pProgMesh->reserveElement(progNumOfElem);//確定した要素数である.

    // --
    // NodeデータをpProgMeshにセット
    // --
    uint i;
    for(i=0; i< mNumOfNode; i++){
        pProgMesh->setNode(mvNode[i]);
    };
    // --
    // Elementデータは,pProgMeshにセットしない.
    //  => 要素については,新たに生成する要素のみがpProgMeshの要素である.
    // --
    

    // ○ ノード(Node)数については,未確定
    //     ->setupEdgeElement,setupFaceElement,setupVolumeNodeでカウントしてその都度セットする.
    // ○ 要素数は確定
    // pProgMesh->setNumOfElement(progNumOfElem);
    //
    //  -> FactoryでのsetupNumOfElement()に変更(vectorのsize()を利用する)


    //debug
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::MWDebug,"prolongation numOfElem",progNumOfElem);
    pLogger->Info(Utility::LoggerMode::MWDebug,"prolongation numOfNode",progNumOfNode);

}


// 各要素の辺に接する要素集合の集計(要素の自分自身のCElement*も集合に含まれる)=> 各要素に辺要素集合をセット
//                && prolongation中間ノードの生成(現LevelのMeshでのノード数に影響なし)
//
// nLevel:2次要素の場合のファーストステップ判定用
//
void CMesh::setupEdgeElement(CMesh *pProgMesh, const uint& nLevel)
{
    // 節点の要素集合は取得済みとする(mvAggregateElement)
    //
    CElement *pElem;
    PairNode edgeNode;
    CNode *pNode0,*pNode1;
    
    // pElemまわりに集められる要素に対する処理用
    //
    CElement *pOtherElem;
    vector<CElement*> vEdgeElem;//一つの要素の辺の要素集合
    uint nOtherEdge;// elem loop外の要素の辺(Edge)index番号

    //vint pairIndex;  pairIndex.resize(2);
    uint numOfEdge;
    uint ielem,iedge;
    uint indexCount;//ノードID生成用
    CNode *inNode;  //中間ノード

    // ID 初期値
    if(pProgMesh){
        //prolongation MeshのノードIDの初期値
        indexCount = pProgMesh->getNodeSize();
    }else{
        //Refine後の2次ノード設定用
        indexCount= mvNode.size();//辺ノードの初期値
    }

    ////debug
    //cout << "indexCount initial = " << indexCount << ", CMesh::setupEdgeElement" << endl;
    
    // 要素 ループ
    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem = mvElement[ielem];
        numOfEdge = pElem->getNumOfEdge();
        
        // 辺ノードを生成してよいか.
        //
        bool bFlag=false;
        if(pElem->getOrder()==ElementOrder::Second && nLevel > 0 ) bFlag=true;
        if(pElem->getOrder()==ElementOrder::First)  bFlag=true;

        // 辺ノードの生成
        if(bFlag){
            // 辺 ループ
            for(iedge=0; iedge< numOfEdge; iedge++){

                edgeNode = pElem->getPairNode(iedge);
                pNode0= edgeNode.first;  pNode1= edgeNode.second;

                // 同一のEdgeを処理することを回避
                if(!pElem->isEdgeElem(pNode0,pNode1)){
                    //pElem->getPairNode(pairIndex,iedge);
                    uint iagg,jagg, elemIndex0, elemIndex1, elemid0, elemid1;
                    uint numOfAggElem0=pNode0->getNumOfAggElem();// ノードの要素集合
                    uint numOfAggElem1=pNode1->getNumOfAggElem();// ノードの要素集合

                    pElem->reserveEdgeElement(iedge,numOfAggElem0);// ノードの要素集合の個数を越えることはないので、ノード要素集合数で確保


                    // 両端のノードの要素集合比較による, 辺の要素集合
                    //
                    for(iagg=0; iagg< numOfAggElem0; iagg++){

                        elemid0= pNode0->getAggElemID(iagg);
                        elemIndex0= moBucket.getIndexElement(elemid0);

                        for(jagg=0; jagg< numOfAggElem1; jagg++){

                            elemid1= pNode1->getAggElemID(jagg);
                            elemIndex1= moBucket.getIndexElement(elemid1);

                            //pElem自身もEdgeElementに入れる.
                            if(elemIndex0 == elemIndex1)
                                pElem->setEdgeElement(iedge, mvElement[elemIndex0]);
                        };
                    };


                    // Edge(辺)中間ノードの生成.
                    inNode= GeneInterNode(pNode0);// 辺ノードの生成

                    // ノード2点間の平均
                    vdouble P0_coord= pNode0->getCoord();
                    vdouble P1_coord= pNode1->getCoord();
                    vdouble In_coord; In_coord.reserve(3);
                    for(uint i=0; i< 3; i++){ In_coord.push_back( (P0_coord[i] + P1_coord[i])*0.5 );}

                    inNode->setCoord(In_coord);//inNode(中間ノード)に座標値をセット
                    inNode->setID(indexCount); //indexCountノード数カウントのセット

                    setupParentNode(pNode0,pNode1,inNode);//prolongater用Nodeのセット

                    pElem->setEdgeInterNode(inNode,iedge);//ループ要素自身に中間ノードをセット
                    pElem->setBoolEdgeElem(pNode0, pNode1);//iedgeにEdgeElementとNodeをセットしたことをスタンプ


                    if(pProgMesh){
                        pProgMesh->setNode(inNode);//上位の頂点ノードに追加
                        //setEdgeNode(inNode); //自分の辺ノードに追加
                    }else{
                        //NULLの場合,Refine後の2次ノード設定
                        //setEdgeNode(inNode);
                    }

                    //-------------------------------------------------
                    // 2次要素の場合：自身のMeshのノード集合に辺ノードをセット
                    //-------------------------------------------------
                    if(pElem->getOrder()==ElementOrder::Second){
                        mvNode.push_back(inNode);
                        mNumOfNode += 1;
                        CAggregateElement *pTAggElem = new CAggregateElement;
                        mvAggElement.push_back( pTAggElem );
                        
                        moBucket.re_resizeBucketNode(indexCount);//新規に追加された節点のID(新MaxID)
                        moBucket.setIndexNode( indexCount, mvNode.size()-1 );
                    }

                    indexCount++;


                    vEdgeElem.clear();
                    // EdgeElement ループ:pElemの辺に集めた要素集合を、要素集合に含まれる要素の各辺にセットする
                    //
                    vEdgeElem= pElem->getEdgeElement(iedge);
                    uint n;
                    for(n=0; n< vEdgeElem.size(); n++){
                        pOtherElem=vEdgeElem[n];

                        //条件1.辺の要素集合にはpElem自身も入っているので,pElemを除外
                        //条件2.辺の要素集合の要素の辺が既にEdge要素をセットしてあるか(isEdgeElem==false)判定
                        if(pOtherElem->getID() != pElem->getID() && !pOtherElem->isEdgeElem(pNode0,pNode1)){
                            nOtherEdge =pOtherElem->getEdgeIndex(pNode0, pNode1);//pNode0,pNode1に該当する辺のIndex番号
                            pOtherElem->setEdgeAggElement(nOtherEdge,vEdgeElem);//集められた要素集合を一度にセット
                            pOtherElem->setBoolEdgeElem(pNode0, pNode1);        //辺の要素集合をセット済みをスタンプ

                            pOtherElem->setEdgeInterNode(inNode,nOtherEdge);//辺集合の要素に中間ノードをセット
                        }
                    };//for(vEdgeElemのループ)

                    //---------------------------------------------------
                    // 2次要素の場合：辺要素集合を"inNode"の要素集合としてセット
                    //---------------------------------------------------
                    if(pElem->getOrder()==ElementOrder::Second && mnSolutionType==SolutionType::FEM){
                        uint index = mvNode.size() - 1;//このルーチンで生成してpushされたばかりなので.size-1
                        CAggregateElement *pAggElem = mvAggElement[index];

                        pAggElem->reserve(vEdgeElem.size());

                        for(uint i=0; i < vEdgeElem.size(); i++){
                            CElement *pEdgeElem = vEdgeElem[i];
                            pAggElem->push(pEdgeElem);
                        };
                    }

                }//ifブロック終端(pElemの一つの辺の処理終端
            };//for ( iedge )
            
        }//if(bGeneEdgeNode) : 辺ノードを生成させるか否か.

    };//for (ielem)

    // 現時点でのノード数をセット
    //pProgMesh->setNumOfNode(indexCount);
    // -> Factoryで行う.

    ////debug
    //cout << "indexCount final   = " << indexCount << ", CMesh::setupEdgeElement" << endl;
    //cout << "pProgMesh NodeSize = " << pProgMesh->getNodeSize() << ", CMesh::setupEdgeElement" << endl;
    
}
//
// 要素が二次要素の場合、要素クラス内で辺ノードをmvNodeに移し替え
//
void CMesh::replaceEdgeNode()
{
    uint ielem;
    // 要素 ループ : 要素内の辺ノードを移し替え
    for(ielem=0; ielem< mNumOfElement; ielem++){
        mvElement[ielem]->replaseEdgeNode();
    };
    
}


// prolongationの為の,中間節点の生成(setupEdgeElementからコール) <= 現LevelのMeshでのノード数に影響なし
//
CNode* CMesh::GeneInterNode(CNode* pNode)
{
    CNode *inNode;
    uint numOfScalar, numOfVector;
    // 要素内では同一種類のノードで構成されるので、片方のノードの型で判断
    // --
    // Scalar
    switch(pNode->getType()){
        case(NodeType::Scalar):
            inNode= new CScalarNode;
            numOfScalar= pNode->numOfScalarParam();
            inNode->resizeScalar(numOfScalar);//DOFの確保

            break;
        case(NodeType::Vector):
            inNode= new CVectorNode;
            numOfVector= pNode->numOfVectorParam();
            inNode->resizeVector(numOfVector);

            break;
        case(NodeType::ScalarVector):
            inNode= new CScalarVectorNode;

            numOfScalar= pNode->numOfScalarParam();
            numOfVector= pNode->numOfVectorParam();

            inNode->resizeScalar(numOfScalar);
            inNode->resizeVector(numOfVector);

            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error,"Node Generation Error, CMesh::GeneInterNode");
            break;
    }

    //inNode->setMGLevel(mMGLevel+1);//prolongation用としてのノードなのでLevelを一段上げる

    ////debug
    //mnDummyCount++;
    //cout << "GeneInterNode count = " << mnDummyCount << endl;

    return inNode;
}

// 各要素の面に隣接する要素を検索 => 要素自身に、面で隣接する要素をセット
//
void CMesh::setupFaceElement(CMesh* pProgMesh)
{
    CElement *pElem;
    vector<CNode*> vFaceCnvNode;
    uint  elemIndex0, elemIndex1, elemIndex2, elemid0, elemid1, elemid2;
    CNode *pNode0,*pNode1,*pNode2;
    CNode *inNode;//中間ノード用ポインター
    uint ielem,isurf,jsurf, iagg,jagg, kagg;
    uint numAggElems0, numAggElems1, numAggElems2;
    uint numOfFace;

    vuint vnShareElems0, vnShareElems1;//共有_要素インデックス

    uint indexCount;// prolongation MeshのノードID
    indexCount = pProgMesh->getNodeSize();


    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem= mvElement[ielem];
        numOfFace= pElem->getNumOfFace();

        for(isurf=0; isurf< numOfFace; isurf++){
            
            vFaceCnvNode = pElem->getFaceCnvNodes(isurf);//面を構成するノードの取得

            pNode0= vFaceCnvNode[0]; pNode1= vFaceCnvNode[1]; pNode2= vFaceCnvNode[2];
            
            //共有 要素インデックス クリア
            vnShareElems0.clear(); vnShareElems1.clear();


            // 既にFaceElemがセットされているか判定(セットされている=> true)
            //
            if(!pElem->isFaceElem(pNode0,pNode1,pNode2)){

                numAggElems0= pNode0->getNumOfAggElem(); vnShareElems0.reserve(numAggElems0);
                numAggElems1= pNode1->getNumOfAggElem();
                //
                // ノード0 とノード1の共有elemIndex => vnShareElems0
                //
                for(iagg=0; iagg< numAggElems0; iagg++){
                for(jagg=0; jagg< numAggElems1; jagg++){
                    elemid0 = pNode0->getAggElemID(iagg);
                    elemIndex0= moBucket.getIndexElement(elemid0);

                    elemid1 = pNode1->getAggElemID(jagg);
                    elemIndex1= moBucket.getIndexElement(elemid1);

                    if(elemIndex0 == elemIndex1)
                        vnShareElems0.push_back(elemIndex0);//ノード0,ノード1 の共有elemIndex

                };//iagg ループ(Vertex周囲の要素インデックス)
                };//jagg ループ

                numAggElems2= pNode2->getNumOfAggElem(); vnShareElems1.reserve(numAggElems2);
                //
                // ノード2 と vnShareElems0(ノード0とノード1の共有要素)との共有elemIndex => vnShareElems1
                //
                for(jagg=0; jagg< vnShareElems0.size(); jagg++){
                for(kagg=0; kagg< numAggElems2; kagg++){
                    elemid2= pNode2->getAggElemID(kagg);
                    elemIndex2= moBucket.getIndexElement(elemid2);
                    if(vnShareElems0[jagg] == elemIndex2)
                        vnShareElems1.push_back(vnShareElems0[jagg]);//ノード0,ノード1,ノード2 の共有elemIndex

                };//kagg ループエンド
                };//jagg ループエンド
                

                // 1.隣接(Adjacent)要素をpElemへセットアップ
                // 2.隣接要素へpElemをセットアップ
                //
                CElement* pAdjElem;

                // 隣接要素が存在する面
                if(vnShareElems1.size()> 1){
                    // pElem自身と隣接要素で2つ. 2つを越える場合は隣接要素形状が2分割されている.<= 現状のメッシュは隣接メッシュは一つとする.
                    for(uint i=0; i< vnShareElems1.size(); i++){
                        
                        uint pElemIndex= moBucket.getIndexElement(pElem->getID());//// '09.10.01

                        if(pElemIndex != vnShareElems1[i]){
                            
                            pAdjElem= mvElement[vnShareElems1[i]];

                            pElem->setFaceElement(pAdjElem, isurf);
                            pElem->setBoolFaceElem(pNode0,pNode1,pNode2);//隣接要素をセットしたことをスタンプ

                            ///////////////////////////
                            // 面ノードの生成,セット
                            inNode= GeneInterNode(pNode0);
                            avgCoord(vFaceCnvNode, inNode);//vFaceCnvNodeの平均座標をセット
                            inNode->setID(indexCount);//IDをセット
                            pProgMesh->setNode(inNode);
                            pElem->setFaceNode(inNode, isurf);// Faceへノードをセット(要素ループしている要素)
                            
                            setupParentNode(vFaceCnvNode,inNode);//prolongater用Nodeのセット
                            
                            //隣接要素の"pNode0,pNode1,pNode2で構成される面へpElemをセット
                            //
                            jsurf= pAdjElem->getFaceIndex(pNode0, pNode1, pNode2);
                            pAdjElem->setFaceElement(pElem, jsurf);
                            pAdjElem->setBoolFaceElem(pNode0,pNode1,pNode2);//隣接要素をセットしたことをスタンプ
                            pAdjElem->setFaceNode(inNode, jsurf);           // Faceへノードをセット(隣の要素)

                            // IDのためのカウントアップ
                            indexCount++;
                        }
                    };
                // 隣接要素が存在しない面
                }else{
                    //////////////////////////////
                    //面ノードの生成,セット
                    inNode= GeneInterNode(pNode0);
                    avgCoord(vFaceCnvNode, inNode);//vFaceCnvNodeの平均座標をセット

                    inNode->setID(indexCount);//IDをセット
                    pProgMesh->setNode(inNode);
                    
                    setupParentNode(vFaceCnvNode,inNode);//prolongater用Nodeのセット

                    // IDのためのカウントアップ
                    indexCount++;
                    
                    pElem->setFaceNode(inNode,isurf);
                    pElem->setBoolFaceElem(pNode0, pNode1, pNode2);
                }
            }//if(!isFaceElem) ブロック エンド
        };//isurf ループ エンド
    };//ielem ループ エンド
    
    // 現時点でのノード数をセット
    // pProgMesh->setNumOfNode(indexCount);
    //   -> FactoryでのsetupNumOfNode()に統一.
}


// 要素の中心にノードを生成、セットアップ
//
void CMesh::setupVolumeNode(CMesh *pProgMesh)
{
    CElement* pElem;
    CNode *pNode,*cntNode;
    vector<CNode*> vLocalNode;
    uint nNumOfVert;
    vdouble vCoord;
    uint ielem, inode;

    uint indexCount;// prolongation MeshのノードID
    indexCount = pProgMesh->getNodeSize();


    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem= mvElement[ielem];

        // Solid要素のみにVolume中心ノードを設置
        //
        //if(pElem->getNumOfFace() > 1){
        if(pElem->getEntityType()==BaseElementType::Solid){
            pNode= pElem->getNode(0);//"局所番号0"のノード <= 生成するノードの種類を決めるために必要.

            cntNode= GeneInterNode(pNode);//要素中心ノードの生成(pNodeと同一種類のノード)
            cntNode->setID(indexCount);

            vCoord.clear();
            vCoord.resize(3);
            vCoord[0]=0.0; vCoord[1]=0.0; vCoord[2]=0.0;
            
            nNumOfVert= pElem->getNumOfVert();// '10.11.24
            vLocalNode.resize(nNumOfVert);

            //局所ノードの座標平均をとる
            for(inode=0; inode< nNumOfVert; inode++){
                pNode= pElem->getNode(inode);
                vLocalNode[inode] = pNode;

                vCoord[0] += pNode->getX();
                vCoord[1] += pNode->getY();
                vCoord[2] += pNode->getZ();
            };
            vCoord[0] /= (double)nNumOfVert;  vCoord[1] /= (double)nNumOfVert;  vCoord[2] /= (double)nNumOfVert;


            //生成した中心ノードに座標値をセット
            cntNode->setCoord(vCoord);
            
            setupParentNode(vLocalNode, cntNode);//prolongater用Nodeのセット

            pProgMesh->setNode(cntNode);
            ++indexCount;

            pElem->setVolumeNode(cntNode);//要素へ中心ノードをセット
        }//if()エンド：要素の面が1を越える=> Solid要素
    };

    // 現時点でのノード数をセット
    // pProgMesh->setNumOfNode(indexCount);
    //   -> FactoryでのsetupNumOfNode()に統一.
}


// 複数のノードの座標平均をpNodeにセット
//
void CMesh::avgCoord(vector<CNode*> vCnvNode, CNode *pNode)
{
    //inNodeの座標を計算(平均値)
    vdouble vCoord; vCoord.resize(3);
    uint i;
    for(i=0; i< 3; i++){ vCoord[i]=0.0;}

    uint numOfFaceNode = vCnvNode.size();
    for(i=0; i< numOfFaceNode; i++){
        vCoord[0] += vCnvNode[i]->getX();
        vCoord[1] += vCnvNode[i]->getY();
        vCoord[2] += vCnvNode[i]->getZ();
    };
    for(i=0; i< 3; i++){
        vCoord[i] /= (double)numOfFaceNode;
    };

    ////debug
    //cout << "vCoord[0]= " << vCoord[0] << ", vCoord[1]= " << vCoord[1] << ", vCoord[2]= " << vCoord[2]<< endl;

    //座標値をinNodeにセットする.
    pNode->setCoord(vCoord);
}


// Refine時に生成される新Node(子Node)の親Nodeを,新Nodeにセットする
// --
// 辺から生成されるNode用
void CMesh::setupParentNode(CNode* pNode0, CNode* pNode1, CNode* inNode)
{
    inNode->reserveParentNode(2);
    inNode->addParentNode(pNode0);
    inNode->addParentNode(pNode1);
}
// 面,要素から生成されるNode用
void CMesh::setupParentNode(vector<CNode*>& vNode, CNode* inNode)
{
    uint numOfParent= vNode.size();
    inNode->reserveParentNode(numOfParent);

    uint ipare;
    for(ipare=0; ipare< numOfParent; ipare++){
        inNode->addParentNode(vNode[ipare]);
    };
}
//// Refine時に生成した子Nodeを,親Nodeにセットする
//// --
//// 辺から生成された子Nodeを両端の親Nodeにセット
//void CMesh::setupChildNode(CNode* pNode0, CNode* pNode1, CNode* inNode)
//{
//    pNode0->addChildNode(inNode);
//    pNode1->addChildNode(inNode);
//}
//// 面(体)から生成された子Nodeを,面(体)を構成する親Nodeにセット (要素)
//void CMesh::setupChildNode(vector<CNode*>& vNode, CNode* inNode)
//{
//    uint numOfParent= vNode.size();
//
//    uint ipare;
//    for(ipare=0; ipare < numOfParent; ipare++){
//        vNode[ipare]->addChildNode(inNode);
//    };
//}


// CommMesh:通信領域
// --
// CommID:通信領域番号
// DommID:計算領域番号
//
// Domm:領域, trasmitDomain:送信先領域
// --
void CMesh::setCommMesh(CCommMesh* pCommMesh)
{
    mvCommMesh.push_back(pCommMesh);

    // CommIDからmvCommMeshのインデックス番号を取得できるようにHashをセット
    uint comID= pCommMesh->getCommID();
    mmCommIndex[comID]= mvCommMesh.size()-1;
}

CCommMesh* CMesh::getCommMesh(const uint& comID)
{
    uint comIndex= mmCommIndex[comID];

    return mvCommMesh[comIndex];
}

// CommMeshのprolongationで不要になったNode,Elementを配列の後ろに移動
//  =>  不要Node=>DNode, 不要Element=>DElement
// --
// 結局このメソッドは,下記の二つの値を取得している;
//
// mNodeEndIndex  <= 計算に使用するNode配列数
// mElemEndIndex  <= 計算に使用するElement配列数
// --
void CMesh::sortMesh()
{
    //debug
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, initial mvNode.size    => ",(uint)mvNode.size());
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, initial mvElement.size => ",(uint)mvElement.size());


    CNode *pNode, *pDNode;
    CElement *pElem, *pDElem;
    uint numOfDNode, numOfDElement;
    
    CCommMesh *pCommMesh;
    uint numOfComm= mvCommMesh.size();
    
    // mvNode,mvElementからDNode,DElementを検索してerase().
    // --
    // ○ mNodeEndIndex <= 計算に使用する配列数 をセット
    // ○ mElemEndIndex <= 計算に使用する配列数 をセット
    
    uint icom;
    uint idel;//eraseしたインデックス管理

    vector<CNode*>::iterator    itNode;
    vector<CElement*>::iterator itElement;

    // 各CommMeshのDNode,DElementを,MeshのmvNode,mvElementから除外.
    //
    for(icom=0; icom< numOfComm; icom++){
        pCommMesh= mvCommMesh[icom];

        numOfDNode= pCommMesh->getNumOfDNode();
        numOfDElement= pCommMesh->getNumOfDCommElement();
      
        //maxIndex= mvNode.size()-1;
        //QuicksortID<CNode*>(mvNode, 0, maxIndex);// クイック・ソート
        sortID<CNode*>(mvNode, mvNode.size());//ソート

        // Nodeの整理(配列のerase,中身(Node*)はそのまま.)
        idel=0;
        for(itNode=mvNode.begin(); itNode< mvNode.end(); itNode++){
            pNode= *itNode;

            // DNodeのerase()
            if(idel < numOfDNode){
                pDNode= pCommMesh->getDNode(idel);

                if(pNode->getID() == pDNode->getID()){
                    mvNode.erase(itNode);
                    idel++; //DNode配列から拾うインデックスを進める.

                    if(itNode != mvNode.begin()) itNode--;
                }
            }
        };//Nodeループ
        
        //maxIndex= CElement.size()-1;
        //QuicksortID<CElement*>(mvElement, 0, maxIndex);// クイック・ソート
        sortID<CElement*>(mvElement, mvElement.size());//ソート

        // Element配列の整理(配列のerase,中身(Element*)はそのまま.)
        idel=0;
        for(itElement=mvElement.begin(); itElement< mvElement.end(); itElement++){
            pElem= *itElement;

            ////debug
            //cout << "pElem ID => " << pElem->getID() << endl;
            
            // DElementのerase()
            if(idel< numOfDElement){
                pDElem= pCommMesh->getDElement(idel);
                if(pDElem->getID() == pElem->getID()){
                    mvElement.erase(itElement);
                    idel++;//DElement配列から拾うインデックスを進める.

                    ////debug
                    //cout << "DElement idel => " << idel << endl;

                    if(itElement != mvElement.begin()) itElement--;
                }
            }
        };//Elementループ
    };//CommMeshループ
    
    // mvNode,mvElementの終端Index(サイズ)を取得
    //  => 計算に使用する配列数
    // --
    mNodeEndIndex= mvNode.size();
    mElemEndIndex= mvElement.size();

    //debug
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mNodeEndIndex => ",mNodeEndIndex);
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mElemEndIndex => ",mElemEndIndex);

    
    // 保全 
    // --
    // mvNode,mvElementにDNode,DElementを追加して全体を保全
    // --
    uint idnode, idelem;
    for(icom=0; icom< numOfComm; icom++){
        pCommMesh= mvCommMesh[icom];
        
        numOfDNode= pCommMesh->getNumOfDNode();
        numOfDElement= pCommMesh->getNumOfDCommElement();
        
        for(idnode=0; idnode< numOfDNode; idnode++){
            pDNode= pCommMesh->getDNode(idnode);
            mvNode.push_back(pDNode);
        };
        for(idelem=0; idelem< numOfDElement; idelem++){
            pDElem= pCommMesh->getDElement(idelem);
            mvElement.push_back(pDElem);
        };
    };

    //debug
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mvNode.size    => ",(uint)mvNode.size());
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mvElement.size => ",(uint)mvElement.size());

    //    --
    //    ↓ 下記処理は,中止 ：階層間でIDが変化すると,不都合
    //    --
    //    // 保全 2
    //    // --
    //    // mvNode,mvElementの配列順序が変わったので,Index番号とIDが不一致
    //    //  => IDとIndex番号が一致するように,ID番号を改変
    //    // --
    //    uint numOfNode= mvNode.size();
    //    uint numOfElement= mvElement.size();
    //    uint index;
    //    for(index=0; index < numOfNode; index++){
    //        pNode= mvNode[index];
    //        pNode->setID(index);
    //    };
    //    for(index=0; index < numOfElement; index++){
    //        pElem= mvElement[index];
    //        pElem->setID(index);
    //    }
}


// CommMesh2 のセット
//
void CMesh::setCommMesh2(CCommMesh2* pCommMesh2)
{
    mvCommMesh2.push_back(pCommMesh2);
    
    mmComm2Index[pCommMesh2->getID()]= mvCommMesh2.size()-1;//Hashデータ:mmComm2Index[commID] => vector Index
}

CCommMesh2* CMesh::getCommMesh2(const uint& comID)
{
    uint index;
    index= mmComm2Index[comID];

    return mvCommMesh2[index];
}



// Refine後の処理
// --
// 辺-面 Node*の削除
// 辺-面 Element*の削除
// --
void CMesh::deleteProgData()
{
    uint ielem;
    for(ielem=0; ielem < mNumOfElement; ielem++)
        mvElement[ielem]->deleteProgData();
}

// vertexのAggElemを削除  <= 全てのMesh処理が終わった後で呼び出す, CMW::FinalizeRefine()
//
void CMesh::deleteAggregate_on_Node()
{
    uint inode;
    for(inode=0; inode < mNumOfNode; inode++)
        mvNode[inode]->deleteAggregate();
}


//
// グループ
//
void CMesh::addElemGrp(CElementGroup* pElemGrp)
{
    mvElementGroup.push_back(pElemGrp);

    uint nGrpID = pElemGrp->getID();
    mmElemGrpID2IX[nGrpID] = mvElementGroup.size() - 1;
}

uint CMesh::getNumOfElemGrp()
{
    return mvElementGroup.size();
}

CElementGroup* CMesh::getElemGrpIX(const uint& index)
{
    return mvElementGroup[index];
}
CElementGroup* CMesh::getElemGrpID(const uint& nGrpID)
{
    uint index = mmElemGrpID2IX[nGrpID];
    return mvElementGroup[index];
}


//////
//////  ID 決定の為のメソッド{ 下位グリッドからの増加節点 }
//////
////uint CMesh::increaseNode(CMesh* pRestMesh)
////{
////    uint nIncreNum = mvNode.size() - pRestMesh->getNumOfNode();
////
////    uint icomm, nNumOfComm = mvCommMesh2.size();
////    CCommMesh2 *pCommMesh, *pRestCommMesh;
////    CHecMPI *pMPI = CHecMPI::Instance();
////
////    uint nNumOfSRComN(0);    //small rank CommNode
////    uint nRestNumOfSRComN(0);//下位(restrict) small rank CommNode
////    uint nIncreComN(0);      //増加したsmall rank CommNode
////    
////    for(icomm=0; icomm < nNumOfComm; icomm++){
////        pCommMesh     = mvCommMesh2[icomm];
////        pRestCommMesh = pRestMesh->getCommMesh2IX(icomm);//下位グリッド通信メッシュ
////
////        // 小Rankとの通信に使用される節点数の増加数
////        if(pCommMesh->getTrasmitRank() <  pMPI->getRank()){
////            nNumOfSRComN = pCommMesh->getCommNodeSize();
////            nRestNumOfSRComN = pRestCommMesh->getCommNodeSize();
////
////            nIncreComN += nNumOfSRComN - nRestNumOfSRComN;
////        }
////    };
////
////    nIncreNum -= nIncreComN;
////
////    return nIncreNum;
////}



