
#include <vector>

//
//  CommHexa.cpp
//
//
//
//                  2009.09.01
//                  2009.09.01
//                  k.Takeda
#include "CommHexa.h"
using namespace pmw;


// construct & destruct
//
CCommHexa::CCommHexa()
{
    //// prolongation-Comm要素 <= rankに所属するか否かは,後に決定 <= 不要な要素も後で利用(REVOCAP対応)
    //mvProgCommElem.reserve(8);

    // Node rank
    mvNodeRank.resize(8);
    mvEdgeRank.resize(12);
    mvFaceRank.resize(6);

    mvbSend.resize(8);
    mvbRecv.resize(8);
    mvbOther.resize(8);

    // Node Index生成したか？Marking
    // DNodeとして取得されたか？Marking
    mvbNodeIXCheck.resize(8);
    mvbDNodeMarking.resize(8);
    uint i;
    for(i=0; i< 8; i++){
        mvbNodeIXCheck[i]=false;
        mvbDNodeMarking[i]=false;

        mvbSend[i]=false;
        mvbRecv[i]=false;
        mvbOther[i]=false;
    }

    // 頂点別の要素集合
    mvvAggCommElem.resize(8);
    mvvNeibCommElemVert.resize(8);
    
    // CommMesh内のグローバルIndex
    mvCommNodeIndex.resize(8);
}

CCommHexa::~CCommHexa()
{
    ;
}


// debug method :所有しているElementの型が一致するか？
//
bool CCommHexa::isTypeCoincidence()
{
    bool bCoin(false);
    
    if(mpElement->getType()==ElementType::Hexa) bCoin=true;
    
    return bCoin;
}


// <<<< CommElmentに移行 >>>>
//
//// prolongation要素の頂点(Node)の計算領域属性
////
//void CCommHexa::setupProgNodeRank(const uint& mgLevel)
//{
//    CEdgeTree *pEdgeTree = CEdgeTree::Instance();// 頂点に接続している辺番号
//    //CFaceTree *pFaceTree = CFaceTree::Instance();// 頂点に接続している面番号
//    CEdgeFaceTree *pEdgeFaceTree = CEdgeFaceTree::Instance();//面を構成する辺番号
//    //CProgElementTree *pProgTree = CProgElementTree::Instance();//progElemの頂点番号の取得
//
//    // 辺,面,体積中心のrank(DomID)を決定
//    // { mvNodeRank <= 頂点ごとのrankはセット済み}
//    // --
//    mvEdgeRank;//辺 ごとの所属計算領域(DomID <=rank)
//    mvFaceRank;//面 ごとの所属計算領域(DomID <=rank)
//    mVolRank;  //体積中心の所属計算領域(DomID <=rank)
//
//    uint iedge,iface;
//    uint *localNodes, *localEdges;
//    int   rankDiff;
//
//    // 辺中心のRank
//    // --
//    for(iedge=0; iedge< 12; iedge++){
//        localNodes= pEdgeTree->getHexaLocalNodeNum(iedge);
//
//        rankDiff= mvNodeRank[localNodes[0]] - mvNodeRank[localNodes[1]];//端点のNodeのrankを比較,小さい方を選択
//        if(mgLevel==1) rankDiff *= -1; //最初のprolongationだけは,逆にRankを決めるための処理
//
//        if(rankDiff <  0) mvEdgeRank[iedge]= mvNodeRank[localNodes[0]];
//        if(rankDiff >= 0) mvEdgeRank[iedge]= mvNodeRank[localNodes[1]];
//    };
//
//    // 面中心のRank
//    // --
//    for(iface=0; iface< 6; iface++){
//        localEdges= pEdgeFaceTree->getHexaFaceConnEdge(iface);
//
//        rankDiff= mvEdgeRank[localEdges[0]] - mvEdgeRank[localEdges[2]];// 対向する辺の値を比較(小さいrankを選択)
//        if(mgLevel==1) rankDiff *= -1;
//
//        if(rankDiff <  0) mvFaceRank[iface]= mvEdgeRank[localEdges[0]];
//        if(rankDiff >= 0) mvFaceRank[iface]= mvEdgeRank[localEdges[2]];
//    };
//
//    // 体積中心のRank
//    // --
//    rankDiff= mvFaceRank[0] - mvFaceRank[1];//対向する面
//    if(mgLevel==1) rankDiff *= -1;
//
//    if(rankDiff < 0) mVolRank= mvFaceRank[0];
//}







//// Send基準のprogElem(prolongationのサブルーチン)
//// --
//void CCommHexa::progSend()
//{
//    CEdgeTree *pEdgeTree = CEdgeTree::Instance();// 頂点に接続している辺番号
//    CFaceTree *pFaceTree = CFaceTree::Instance();// 頂点に接続している面番号
//    CProgElementTree *pProgTree = CProgElementTree::Instance();//progElemの頂点番号の取得
//
//    uint ivert,iedge,iface;
//    uint *pconnEdge, *pconnFace;
//    CCommElement* pProgCommElem;
//    CElement* pProgElem;
//    uint progvert;//progElemの頂点番号
//
//    JudgeSend();//中間ノードのSend判定
//
//    // ProgCommElementの生成
//    // --
//    // 頂点が"Send属性"の場所に作られる分割要素を,ProgCommElemとする
//    for(ivert=0; ivert< 8; ivert++){
//        if(mvbSendVertNode[ivert]){
//            pProgElem= mpElement->getProgElem(ivert);// 頂点番号順にProgElemが入っている.
//
//            //-- CommHexa生成
//            pProgCommElem= new CCommHexa();
//            pProgCommElem->setElement(pProgElem);
//            mvProgCommElem.push_back(pProgCommElem);
//            //--
//            pconnEdge= pEdgeTree->getHexaConnEdge(ivert);// 頂点に接続している辺番号
//            pconnFace= pFaceTree->getHexaConnFace(ivert);// 頂点に接続している面番号
//
//            // ProgCommElemの頂点のSend-Recv属性
//            // --
//            pProgCommElem->setSendVert(ivert);
//
//            for(iedge=0; iedge< 3;iedge++){
//                if(mvbSendEdgeNode[pconnEdge[iedge]]){
//                    progvert= pProgTree->getHexaEdgeProgVert(pconnEdge[iedge],ivert);
//                    pProgCommElem->setSendVert(progvert);
//                }
//            }
//            for(iface=0; iface< 3; iface++){
//                if(mvbSendFaceNode[pconnFace[iface]]){
//                    progvert= pProgTree->getHexaFaceProgVert(pconnFace[iface],ivert);
//                    pProgCommElem->setSendVert(progvert);
//                }
//            }
//            progvert= pProgTree->getHexaVolProgVert(ivert);
//            pProgCommElem->setRecvVert(progvert);//頂点の属性の反対の属性
//        }
//    };
//}
//
//// Recv基準のprogElem(prolongationのサブルーチン)
//// --
//void CCommHexa::progRecv()
//{
//    CEdgeTree *pEdgeTree = CEdgeTree::Instance();// 頂点に接続している辺番号
//    CFaceTree *pFaceTree = CFaceTree::Instance();// 頂点に接続している面番号
//    CProgElementTree *pProgTree = CProgElementTree::Instance();//progElemの頂点番号の取得
//
//    uint ivert,iedge,iface;
//    uint *pconnEdge, *pconnFace;
//    CCommElement* pProgCommElem;
//    CElement* pProgElem;
//    uint progvert;//progElemの頂点番号
//
//
//    JudgeRecv();//中間ノードのRecv判定
//
//
//    // ProgCommElemの生成
//    // --
//    // 頂点が"Recv属性"の場所に作られる分割要素を,ProgCommElemとする
//    for(ivert=0; ivert< 8; ivert++){
//        if(mvbRecvVertNode[ivert]){
//            pProgElem= mpElement->getProgElem(ivert);//Hexaは,頂点番号順にProgElemが入っている.
//
//            //-- CommHexa生成
//            pProgCommElem= new CCommHexa();
//            pProgCommElem->setElement(pProgElem);
//            mvProgCommElem.push_back(pProgCommElem);
//            //--
//            pconnEdge= pEdgeTree->getHexaConnEdge(ivert);// 頂点に接続している辺番号
//            pconnFace= pFaceTree->getHexaConnFace(ivert);// 頂点に接続している面番号
//
//            // ProgCommElemの頂点のSend-Recv属性
//            // --
//            pProgCommElem->setRecvVert(ivert);
//
//            for(iedge=0; iedge< 3;iedge++){
//                if(mvbRecvEdgeNode[pconnEdge[iedge]]){
//                    progvert= pProgTree->getHexaEdgeProgVert(pconnEdge[iedge],ivert);
//                    pProgCommElem->setRecvVert(progvert);
//                }
//            }
//            for(iface=0; iface< 3; iface++){
//                if(mvbRecvFaceNode[pconnFace[iface]]){
//                    progvert= pProgTree->getHexaFaceProgVert(pconnFace[iface],ivert);
//                    pProgCommElem->setRecvVert(progvert);
//                }
//            }
//            progvert= pProgTree->getHexaVolProgVert(ivert);
//            pProgCommElem->setSendVert(progvert);//頂点の属性の反対の属性
//        }
//    };
//}
















