//
//  CommunicationTable.cpp
//
//
//
//                      2009.06.19
//                      2009.06.19
//                      k.Takeda
#include <vector>

#include "Element.h"

#include "CommElement.h"
#include "CommunicationMesh.h"
using namespace pmw;

#include <iostream>
// construct & destruct
//
CCommMesh::CCommMesh(CIndexBucket* pBucket)
{
    mpBucket = pBucket;
    mpLogger= Utility::CLogger::Instance();
}
CCommMesh::~CCommMesh()
{
    //debug
    std::cout << "~CCommMesh" << std::endl;
}

// メッシュの再分割:prolongation(refineMesh)
// --
// ToDo:
// --
// mvCommElementの収集
// mvNodeの収集(CommMesh内のIndex番号)
// (CommMeshのIndex)と(MeshのIndex)のmap
// mvSendNode,mvRecvNodeの収集 => MeshのIndexとのmap
// mvDNode,mvDElementの収集

//// prolongation のための辺,面,体積中心のNodeのランク設定
//// => Factoryから直接CommElementのsetupProgNodeRank()が呼ばれているので,CommMeshには無い.



// CommElementの配分(Allocate)
//  => 通信:CommElement,非通信:DCommElement,全体:CommElementAll
// --
//  prolongation時(Refine)の,Communication(CommMesh)内でのIndex生成(Global番号相当)
//  =>
//   CommElementID, CommNodeIDの生成(全てのCommElementに連番：CommとDCommを分けずにIDを打つ)
//   ElementにCommElementのIDをセット, ElementにCommElementに含まれることをスタンプ
// --
void CCommMesh::AllocateCommElement()
{
    uint numOfCommElemAll= mvCommElementAll.size();
    CCommElement* pCommElem;
    CElement* pElem;
    uint icom;
    for(icom=0; icom< numOfCommElemAll; icom++){
        pCommElem= mvCommElementAll[icom];
        
        pCommElem->sortNodeRank(mRankID, mTransmitRankID);

        pElem= pCommElem->getElement();

        if(pCommElem->isCommElement()){

            mvCommElement.push_back(pCommElem);
            pElem->interCommElem(); //ElementにCommElementに含まれている事をスタンプ

        }else if(!pCommElem->isRCommElement()){//通信に使われず&& 計算領域に含まれない(RCommElement は,頂点が全てmyRank)

            mvDCommElement.push_back(pCommElem);
            pElem->interDCommElem();//ElementにDCommElementに含まれていることをスタンプ
        }else{
            pElem->interRCommElem();//CommMeshに含まれるが,通信には使わず(全てがmyRank),DCommElementでもない.
        }

        pCommElem->setID(icom);
        pElem->setCommID(icom);

        ////debug
        //cout << "pElem ID => " << pElem->getID() << ",AllocateCommElement CommElementID =>" << icom << endl;
    };
}

// 引数:vElementは,Mesh全体のElements
//
void CCommMesh::setupAggCommElement(vector<CElement*> &vElement)
{
    uint numOfCommElem= mvCommElementAll.size();
    uint numOfVert;
    uint numOfAggElem;
    CElement* pElem;
    CCommElement* pCommElem;
    CCommElement* pNeibCommElem;
    uint neibComID, neibVert;
    CNode* pNode;
    uint icome,ivert,iagg;
    uint elemID,elemIndex;

    for(icome=0; icome< numOfCommElem; icome++){
        pCommElem= mvCommElementAll[icome];

        numOfVert= pCommElem->getNumOfVert();
        //--
        for(ivert=0; ivert< numOfVert; ivert++){
            pNode= pCommElem->getNode(ivert);

            numOfAggElem= pNode->getNumOfAggElem();
            
            //ElemからNodeの周辺の要素を取得 => CommElemにセット
            //--
            for(iagg=0; iagg< numOfAggElem; iagg++){
                
                elemID= pNode->getAggElemID(iagg);//頂点で接続している要素ID
                elemIndex= mpBucket->getIndexElement(elemID);

                pElem= vElement[elemIndex];
                
                // 取得したElementが,CommMeshに該当すれば,要素を集める.
                if(pElem->isInterCommElem()||pElem->isInterDCommElem()||pElem->isInterRCommElem()){

                    neibComID= pElem->getCommID();

                    //頂点で接続しているCommElem
                    pNeibCommElem= mvCommElementAll[neibComID];
                    pCommElem->setAggCommElement(ivert, pNeibCommElem);

                    //接続先の頂点番号
                    neibVert= pNode->getNeibElemIDVert(elemID);
                    pCommElem->setNeibCommElemVert(ivert,neibVert);
                }
            };//iagg ループ
        };//ivert ループ
    };//icome ループ
}

// CommElementAllから,mvNode全体,SendNode,RecvNodeを取得
// --
// => Send,Recv収集
//   DCommElement,DNode のソート
//   => Mesh本体のElement,Nodeの配列ソートに使用(使用するMeshと,カップラー向けのみのMeshに振り分け)
// --
void CCommMesh::sortCommNodeIndex()
{
    uint numOfCommElem= mvCommElementAll.size();
    uint numOfVert;
    CCommElement* pCommElem;
    uint icome, ivert;

    uint comNodeIndex(0);//CommNodeIndexカウンター

    // mvCommElementAllからNode全体を取得 => 全体にCommNodeIndexをつけていく. => Nodeのランクも収集
    // --
    for(icome=0; icome< numOfCommElem; icome++){
        pCommElem= mvCommElementAll[icome];
        
        numOfVert= pCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++){
            // CommElement内部で重複しないようにマーキングしながらcomNodeIndexをカウントアップ
            //   => Node*をmvNodeに収集. Nodeランクも収集
            //
            pCommElem->setCommNodeIndex(ivert, comNodeIndex, mvNode, mvNodeRank);
        };
    };
    //debug
    mpLogger->Info(Utility::LoggerMode::Debug, "CommMesh::sortCommNodeIndex mvNode.size => ", (uint)mvNode.size());
    
    
    
    // Send,Recvノードの収集
    // --
    // Indexが重複しないように,スタンプを用意
    // --
    vbool vIndexCheck;  uint numOfNode= mvNode.size();
    vIndexCheck.reserve(numOfNode);
    for(uint i=0; i< numOfNode; i++){
        vIndexCheck.push_back(false);
    }
    uint nIndex, rank;
    numOfCommElem= mvCommElement.size();

    for(icome=0; icome< numOfCommElem; icome++){
        pCommElem= mvCommElement[icome];
        
        numOfVert= pCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++){

            nIndex= pCommElem->getCommNodeIndex(ivert);//CommElementの頂点のCommNodeIDを取得

            if(!vIndexCheck[nIndex]){
                rank= mvNodeRank[nIndex];

                if(mRankID==rank){
                    mvSendNode.push_back(mvNode[nIndex]);
                    mvSendCommNodeID.push_back(nIndex);
                }
                if(mTransmitRankID==rank){
                    mvRecvNode.push_back(mvNode[nIndex]);
                    mvRecvCommNodeID.push_back(nIndex);
                }
                vIndexCheck[nIndex]= true;//取得済み
            }
        };
    };

    
    // DCommElementからElementを取得して,Mesh内の計算不要Elementを取得
    // 計算不要Elementをソート => MeshでのElement並び替えに使用
    uint numOfDCommElem= mvDCommElement.size();
    CCommElement* pDCommElem;
    CElement* pDElem;
    for(icome=0; icome< numOfDCommElem; icome++){
        pDCommElem= mvDCommElement[icome];
        pDElem = pDCommElem->getElement();

        mvDElement.push_back(pDElem);
    };
    //debug
    mpLogger->Info(Utility::LoggerMode::Debug, "CommMesh::sortCommNodeIndex mvDElement.size => ", (uint)mvDElement.size());
    
    // DCommElementからDNodeを,重複なきようにスタンプを押しながら取得
    //
    for(icome=0; icome< numOfDCommElem; icome++){
        pDCommElem= mvDCommElement[icome];
        
        numOfVert= pDCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++){
            //// ↓ ////
            pDCommElem->getDNode(ivert, mvDNode);//内部でDNodeか否かの判定処理を行って,mvDNodeにセットする.
        };
    };
    //debug
    mpLogger->Info(Utility::LoggerMode::Debug, "CommMesh::sortCommNodeIndex mvDNode.size => ", (uint)mvDNode.size());


    // 計算不要Node, Elementのソート
    // --
    // mvDNodeをソート => MeshでのNode並び替えに使用
    // mvDElementをソート => MeshでのElement並び替えに使用
    // ----
    uint maxIndex;
    
    if(mvDNode.size() > 0){
        //maxIndex= mvDNode.size()-1;
        //QuicksortID<CNode*>(mvDNode, 0, maxIndex);// クイック・ソート
        sortID<CNode*>(mvDNode, mvDNode.size());//ソート
    }
    if(mvDElement.size() > 0){
        //maxIndex= mvDElement.size()-1;
        //QuicksortID<CElement*>(mvDElement, 0, maxIndex);// クイック・ソート
        sortID<CElement*>(mvDElement, mvDElement.size());// ソート
    }
    ////debug
    //for(uint i=0; i< mvDElement.size(); i++){
    //    mpLogger->Info(Utility::LoggerMode::Debug,"CommMesh::sortCommNodeIndex mvDElement ID => ", (uint)mvDElement[i]->getID());
    //};
}

// 
// mapデータのセットアップ
//  mmCommNodeIX : NodeID => CommNodeID(Index番号)のHash
//  mmCommElementIX : ElementID => CommElementID(Index番号)のHash
// --
void CCommMesh::setupMapID2CommID()
{
    uint numOfCommElem= mvCommElementAll.size();
    uint numOfCommNode= mvNode.size();

    CNode* pNode;
    CElement* pElem;
    CCommElement* pCommElem;

    uint index;
    // mapデータのセット
    // --
    // Node-ID から CommMeshのNodeインデックス
    // --
    for(index=0; index< numOfCommNode; index++){
        pNode= mvNode[index];

        mmCommNodeIX[pNode->getID()]= index;
    };
    // Element-ID から CommMeshのCommElementインデックス
    // --
    for(index=0; index< numOfCommElem; index++){
        pCommElem= mvCommElementAll[index];
        pElem= pCommElem->getElement();

        mmCommElementIX[pElem->getID()]= index;
    };
}




























