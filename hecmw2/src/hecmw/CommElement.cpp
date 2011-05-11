
#include <vector>

#include "Element.h"

//
//  CommElement.cpp
//
//  通信テーブルオブジェクト(CommTable)用途
//  通信インデックス管理用 要素
//
//
//                  2009.08.28
//                  2009.08.28
//                  k.Takeda
#include "CommElement.h"
using namespace pmw;


// construct & destruct
//
CCommElement::CCommElement()
{
    mpLogger = Utility::CLogger::Instance();
    
    mbCommunication = false;
    mbRCommunication= false;
}

CCommElement::~CCommElement()
{
    // Send-Recv 局所Node番号
    // --
    //
    delete[] mvbSend; //頂点番号のノードがSend？ => CommMeshのSendノード収集に使用
    delete[] mvbRecv; //頂点番号のノードがRecv？ => CommMeshのRecvノード収集に使用
    delete[] mvbOther;//頂点番号のノードは無関係なRank？

    // Node選択されたかマーキングするためのboolean
    // --
    delete[] mvbNodeIXCheck; //グローバルIndex生成時に,Nodeが選択されたかマーキング
    delete[] mvbDNodeMarking;//DNodeとして選ばれたかマーキング
}


// prolongation要素の頂点(Node)の計算領域属性
// --
// 辺,面,体積中心の rank を決定{ mvNodeRank <= 頂点ごとのrankはセット済み}
// --
void CCommElement::setupProgNodeRank(const uiint& mgLevel)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();// 頂点に接続している辺番号
    //CFaceTree *pFaceTree = CFaceTree::Instance();// 頂点に接続している面番号
    CEdgeFaceTree *pEdgeFaceTree = CEdgeFaceTree::Instance();//面を構成する辺番号
    //CProgElementTree *pProgTree = CProgElementTree::Instance();//progElemの頂点番号の取得

    
    // コンストラクタでサイズは決定されている.
    // --
    uiint numOfEdge= mvEdgeRank.size();//辺 ごとの所属計算領域(rank)
    uiint numOfFace= mvFaceRank.size();//面 ごとの所属計算領域(rank)
    
    
    uiint iedge,iface;
    uiint *localNodes, *localEdges;
    int   rankDiff;

    // 辺中心のRank
    // --
    for(iedge=0; iedge< numOfEdge; iedge++){
        localNodes= pEdgeTree->getHexaLocalNodeNum(iedge);
        
        rankDiff= mvNodeRank[localNodes[0]] - mvNodeRank[localNodes[1]];//端点のNodeのrankを比較,小さい方を選択
        if(mgLevel==1) rankDiff *= -1; //最初のprolongationだけは,逆にRankを決めるための処理
        
        if(rankDiff <  0) mvEdgeRank[iedge]= mvNodeRank[localNodes[0]];
        if(rankDiff >= 0) mvEdgeRank[iedge]= mvNodeRank[localNodes[1]];
    };
    ////debug
    //for(iedge=0; iedge< numOfEdge; iedge++)
    //    mpLogger->Info(Utility::LoggerMode::Debug, "CommElement::setupProgNodeRank, mvEdgeRank=>", (uint)mvEdgeRank[iedge]);

    
    // 面中心のRank
    // --
    for(iface=0; iface< numOfFace; iface++){
        localEdges= pEdgeFaceTree->getHexaFaceConnEdge(iface);
        
        rankDiff= mvEdgeRank[localEdges[0]] - mvEdgeRank[localEdges[2]];// 対向する辺の値を比較(小さいrankを選択)
        if(mgLevel==1) rankDiff *= -1;
        
        if(rankDiff <  0) mvFaceRank[iface]= mvEdgeRank[localEdges[0]];
        if(rankDiff >= 0) mvFaceRank[iface]= mvEdgeRank[localEdges[2]];
    };
    ////debug
    //for(iface=0; iface< numOfFace; iface++)
    //    mpLogger->Info(Utility::LoggerMode::Debug, "CommElement::setupProgNodeRank, mvFaceRank=>", (uint)mvFaceRank[iface]);


    // 体積中心のRank
    // --
    if(numOfFace > 4){
        rankDiff= mvFaceRank[0] - mvFaceRank[1];//対向する面
        if(mgLevel==1) rankDiff *= -1;

        if(rankDiff < 0) mVolRank= mvFaceRank[0];
        if(rankDiff >=0) mVolRank= mvFaceRank[1];

    }else if(numOfFace == 1){
        mVolRank= mvFaceRank[0];//Quad,Triangle
    }else{
        mVolRank= mvEdgeRank[0];//Beam
    }
    ////debug
    //mpLogger->Info(Utility::LoggerMode::Debug,"CommElement::setupProgNodeRank, mVolRank=>",(uint)mVolRank);

}


// 通信に用いる要素か否か
// --
void CCommElement::sortNodeRank(const uiint& myRank, const uiint& transRank)
{
    //int  rankDiff;
    //int  justRank = transRank - myRank;
    mvSendNode.clear(); mvRecvNode.clear(); mvOtherNode.clear();

    CNode* pNode;
    uiint numOfVert(mpElement->getNumOfNode());
    uiint ivert;

    for(ivert=0; ivert< numOfVert; ivert++){
        pNode= mpElement->getNode(ivert);
        
        if(mvNodeRank[ivert]==myRank){
            mvbSend[ivert]=true;// CommMeshのSendNode収集に使用
            mvSendNode.push_back(pNode);
        }
        else if(mvNodeRank[ivert]==transRank){
            mvbRecv[ivert]=true;// CommMeshのRecvNode収集に使用
            mvRecvNode.push_back(pNode);
        }
        else{
            mvbOther[ivert]=true;
            mvOtherNode.push_back(pNode);
        }
    };
    // 通信に使用されるCommElementか否かの判定.
    // --
    if(mvOtherNode.size()==0 && mvSendNode.size() > 0 && mvRecvNode.size() > 0){
        mbCommunication= true;
    }
    // 全ての頂点が送信(myRank)の場合は,通信に用いられないが,通常の計算領域に含まれる
    //  => RCommElement(CommMeshには含まれ:1.通信には用いられない. 2.計算には用いる)
    // --
    if(mvSendNode.size() == numOfVert){
        mbRCommunication= true;
    }

//    //debug
//    cout << "CommElement sortNodeRank, mpElement->ID " << mpElement->getID()
//         << ", mbCommunication:"  << mbCommunication
//         << ", mbRCommunication:" << mbRCommunication
//         << endl;
}

// 隣接CommElemへのCommNodeIndex割り振り
//  => setCommNodeIndexからコール
//
void CCommElement::setNeibCommNodeIndex(const uiint& ivert, const uiint& comID)
{
    if(!mvbNodeIXCheck[ivert]){
        mvCommNodeIndex[ivert]= comID;
        mvbNodeIXCheck[ivert]= true;
    }
}

// CommMesh内のグローバルNodeインデックス(comIDをメソッド内でカウントアップ)
//  -> CommMeshの"sortCommNodeIndex"のループにインデックスの初期値が存在している.
//  -> 引数:vNode == CommMeshのmvNode
// --
void CCommElement::setCommNodeIndex(const uiint& ivert, uiint& comID, vector<CNode*>& vNode, vuint& vCommMeshNodeRank)
{
    CNode* pNode;
    // markingされていないノードにIndex割り振り,Indexのカウントアップ
    //--
    if(!mvbNodeIXCheck[ivert]){

        ////debug
        //cout << "CommElement::setCommNodeIndex comID => " << comID << endl;
        
        mvCommNodeIndex[ivert]= comID;
        mvbNodeIXCheck[ivert] = true;// marking (スタンプ押し)
        
        // CommMesh本体へ選択したNodeをpush.
        // --
        pNode = mpElement->getNode(ivert);
        vNode.push_back(pNode);//<<<<<<<<-- CommMeshのmvNode

        // CommMesh本体のNodeRank配列を生成
        // --
        vCommMeshNodeRank.push_back(mvNodeRank[ivert]);

        
        // 頂点接続しているCommElement関係
        //
        vector<CCommElement*> vNeibCommElem= mvvAggCommElem[ivert];
        vuint vNeibVert = mvvNeibCommElemVert[ivert];
        
        CCommElement* neibCommElem;
        uiint neibVert;
        
        uiint numOfAgg= vNeibCommElem.size();
        uiint iagg;
        // 頂点を共有しているCommElementのCommNodeにもマーキング&&Indexセット
        // --
        for(iagg=0; iagg< numOfAgg; iagg++){
            neibCommElem= vNeibCommElem[iagg];
            neibVert= vNeibVert[iagg];

            neibCommElem->setNeibCommNodeIndex(neibVert, comID);//頂点接続しているCommElemにもComNodeIDをセット
        };
        // インデックスのカウントアップ
        comID += 1;
    }
}


// DNode判定後に提供
// --
// => mbCommunication==false  但し-> Nodeを共有しているCommElemが,mbCommunication==trueの場合は除外
//
void CCommElement::getDNode(const uiint& ivert, vector<CNode*>& vDNode)
{
    CNode* pDNode;

    if(!mbCommunication){
    if(!mvbDNodeMarking[ivert]){
        
        // 頂点接続しているCommElement関係
        vector<CCommElement*> vNeibCommElem= mvvAggCommElem[ivert];
        vuint vNeibVert = mvvNeibCommElemVert[ivert];
        
        CCommElement* neibCommElem;
        uiint neibVert;

        ////debug
        //mpLogger->Info(Utility::LoggerMode::Debug,"CommElement::getDNode, vNeibCommElem.size => ",(uint)vNeibCommElem.size());
        //mpLogger->Info(Utility::LoggerMode::Debug,"CommElement::getDNode, vNeibVert.size => ",(uint)vNeibVert.size());

        uiint numOfAgg= vNeibCommElem.size();
        uiint iagg;
        bool bCommElem(false);
        // 頂点を共有しているCommElementの中に通信に使用されるCommElementがあるか?
        //   => isCommElementが一つでもtrueならば,通信に使用される.
        // --
        for(iagg=0; iagg< numOfAgg; iagg++){
            neibCommElem= vNeibCommElem[iagg];
            
            if(neibCommElem->isCommElement()){
                bCommElem=true;
            }
        };//頂点集合要素のループ


        // 頂点の周囲のCommElementの中に,通信CommElementは存在しない.
        //  => DNodeである.
        if(!bCommElem){
            pDNode = mpElement->getNode(ivert);
            markingDNode(ivert);//mvbDNodeMarking[ivert]=true;

            //////引数のvDNodeにセットアップ
            vDNode.push_back(pDNode);

            // 頂点を共有しているCommElementのNodeにもマーキング
            // --
            for(iagg=0; iagg< numOfAgg; iagg++){
                neibCommElem= vNeibCommElem[iagg];
                neibVert= vNeibVert[iagg];

                if(!neibCommElem->isMarkingDNode(neibVert)){
                    neibCommElem->markingDNode(neibVert);
                }
            };//頂点集合要素のループ
        }
        
    }//if(!mvbDNodeMarking)
    }//if(!mbCommunication)
}

// DNodeとしてカウント済みであることをマーキング
//
void CCommElement::markingDNode(const uiint& ivert)
{
    mvbDNodeMarking[ivert]=true;
}



















