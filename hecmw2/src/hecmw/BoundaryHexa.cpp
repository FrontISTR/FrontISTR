//
// BoundaryHexa.cpp
//
//          2010.07.02
//          k.Takeda
#include <vector>

#include "BoundaryHexa.h"
#include "ElementProperty.h"
using namespace pmw;


//uint CBoundaryHexa::mnElemType = ElementType::Hexa;
//uint CBoundaryHexa::mNumOfFace = NumberOfFace::Hexa();
//uint CBoundaryHexa::mNumOfEdge = NumberOfEdge::Hexa();
//uint CBoundaryHexa::mNumOfNode = NumberOfVertex::Hexa();

CBoundaryHexa::CBoundaryHexa()
{
    // Volume形状
    // ----
    // mvbMarkingEdge初期化
    // mvbMarkingFace初期化
    // ----
    uiint i;
    // 辺
    mvbMarkingEdge = new bool[NumberOfEdge::Hexa()];
    for(i=0; i < NumberOfEdge::Hexa(); i++){
        mvbMarkingEdge[i]=false;
    };
    mvEdgeBNode.resize(NumberOfEdge::Hexa());
    mvEdgeNeibVol.resize(NumberOfEdge::Hexa());

    // 面
    mvbMarkingFace = new bool[NumberOfFace::Hexa()];
    for(i=0; i < NumberOfFace::Hexa(); i++){
        mvbMarkingFace[i]=false;
    };
    mvFaceBNode.resize(NumberOfFace::Hexa());
    mvFaceNeibVol.resize(NumberOfFace::Hexa());

    for(i=0; i < NumberOfFace::Hexa(); i++){
        mvFaceBNode[i]= NULL;
    }
    mpVolBNode= NULL;
}
CBoundaryHexa::~CBoundaryHexa()
{
    ;
}


// Volume 形状情報
//
uiint CBoundaryHexa::getElemType()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    
    switch(mnOrder){
        case(ElementOrder::First):
            return ElementType::Hexa;
        case(ElementOrder::Second):
            return ElementType::Hexa2;
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryHexa::getElemType, invalid mnOrder");
            return ElementType::Limit;
    }
}
uiint CBoundaryHexa::getNumOfEdge()
{
    return NumberOfEdge::Hexa();
}
uiint CBoundaryHexa::getNumOfFace()
{
    return NumberOfFace::Hexa();
}
uiint CBoundaryHexa::getNumOfNode()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();

    switch(mnOrder){
        case(ElementOrder::First):
            return NumberOfNode::Hexa();
        case(ElementOrder::Second):
            return NumberOfNode::Hexa2();
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryHexa::getNumOfNode, invalid mnOrder");
            return NumberOfNode::Default();
    }
}
uiint CBoundaryHexa::getNumOfVert()
{
    return NumberOfVertex::Hexa();
}

void CBoundaryHexa::setOrder(const uiint& order)
{
    mnOrder = order;
    
    switch(mnOrder){
    case(ElementOrder::First):
        mvBNode.resize(NumberOfNode::Hexa());
        break;

    case(ElementOrder::Second):
        mvBNode.resize(NumberOfNode::Hexa2());
        break;
    }
}

// Edge番号 -> pariBNode
//
PairBNode CBoundaryHexa::getPairBNode(const uiint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uiint *pLocalNum;
    
    pLocalNum= pEdgeTree->getHexaLocalNodeNum(iedge);

    uiint index1st = pLocalNum[0];
    uiint index2nd = pLocalNum[1];

    PairBNode pairBNode;
    pairBNode.first = mvBNode[index1st];
    pairBNode.second= mvBNode[index2nd];

////    // debug
////    CNode *pNode0, *pNode1;
////    pNode0 = mvBNode[index1st]->getNode(); pNode1 = mvBNode[index2nd]->getNode();
////    uint nElemEdge = mpElement->getEdgeIndex(pNode0, pNode1);
////    cout << "BoundaryHexa::getPairBNode, Mesh-Elem  Edge=" << nElemEdge << ", iedge=" << iedge << endl;

    return pairBNode;
}


// pairBNode -> Edge番号
//
uiint& CBoundaryHexa::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uiint ivert= getVertIndex(pairBNode.first);
    uiint jvert= getVertIndex(pairBNode.second);
    
    //cout << "BoundaryHexa::getEdgeID, ivert " << ivert << ", jvert " << jvert << endl;

    return pEdgeTree->getHexaEdgeIndex(ivert,jvert);
}


// 面を構成するBNode配列
//
vector<CBoundaryNode*> CBoundaryHexa::getFaceCnvNodes(const uiint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    CBoundaryNode *pBNode;
    uiint  ivert, i, numOfVert;
    uiint *pvIndex;
    
    pvIndex= pFaceTree->getLocalNodeHexaFace(iface);
    numOfVert= pFaceTree->getHexaFaceNumOfVert(iface);

    vector<CBoundaryNode*> vFaceCnvNodes;//面を構成するBNode
    vFaceCnvNodes.reserve(numOfVert);

    for(i=0; i < numOfVert; i++){
        ivert= pvIndex[i];
        pBNode= mvBNode[ivert];
        vFaceCnvNodes.push_back(pBNode);
    }

    return vFaceCnvNodes;
}

// BNode配列 -> Face番号
//
uiint& CBoundaryHexa::getFaceID(vector<CBoundaryNode*>& vBNode)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();

    vuint vLocalVert;
    CBoundaryNode *pBNode;
    uiint ibnode, numOfBNode= vBNode.size();

    vLocalVert.reserve(4);
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= vBNode[ibnode];
        vLocalVert.push_back(getVertIndex(pBNode));
    };

    return pFaceTree->getHexaFaceIndex2(vLocalVert);
}

// 辺を構成する頂点番号
uiint* CBoundaryHexa::getLocalNode_Edge(const uiint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    return pEdgeTree->getHexaLocalNodeNum(iedge);
}
// 面を構成する頂点番号
uiint* CBoundaryHexa::getLocalNode_Face(const uiint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();

    return pFaceTree->getLocalNodeHexaFace(iface);
}



// Refine Volume再分割
// ----
void CBoundaryHexa::refine(uiint& countID, const vuint& vDOF)
{
    CBoundaryVolume *pProgVol;
    CElement   *pProgElem;
    CNode      *pNode;
    CBoundaryNode *pBNode;
    uiint iprog, numOfProg;
    uiint nProgPos;
    double progCubicVol;//体積
    double coef;

    numOfProg= 8;
    mvProgVolume.reserve(numOfProg);

    for(iprog=0; iprog < numOfProg; iprog++){

        pProgVol = new CBoundaryHexa;// <<<<<<<<<<<<<<<< new

        pProgVol->setID(countID);
        pProgVol->setOrder(mnOrder);//resizeBNodeも含む
        countID++;

        mvProgVolume.push_back(pProgVol);

        nProgPos = dividHexa(iprog, pProgVol);//子供がぶら下がっている場所::nProgPosはBNodeの頂点番号(配列番号)

        //子要素を取得
        pBNode= mvBNode[nProgPos];
        pNode = pBNode->getNode();
        //
        //// ↓ これでもOK
        ////uint eVert = mpElement->getLocalVertNum(pNode->getID());
        ////pProgElem= mpElement->getProgElem(eVert);// <<<<<<<<<<<<<< Nodeの頂点番号にぶらさがっている
        pProgElem = mpElement->getProgElem_NodeID(pNode->getID());//<<<<<<<<<<<<<< Nodeにぶらさがっている(Element内部で頂点に変換している)

        // 子Volに子要素をセット
        pProgVol->setElement(pProgElem);
        pProgVol->setElementID(pProgElem->getID());


        // 体積比による境界値の分配
        // ----
        progCubicVol= pProgVol->calcVolume();
        coef = progCubicVol/mCubicVolume;

        ////debug
        //cout << "BoundaryHexa::refine, progCubicVol==" << progCubicVol << ", mCubicVol==" << mCubicVolume << endl;

        distValue(pProgVol, coef, vDOF);

    };// iprog loop_end
        
}

// BoundaryVolumeの体積
// ----
double& CBoundaryHexa::calcVolume()
{
    CDiscreteVolume *pDiscrete= CDiscreteVolume::Instance();
    uiint* discreHexa;
    uiint  i,ii;
    CNode* vNode[4];

    mCubicVolume= 0.0;

    //6個のTetraの体積を加算
    for(i=0; i< 6; i++){
        discreHexa= pDiscrete->getHexaDiscrete(i);
        for(ii=0; ii< 4; ii++){
            vNode[ii]= mvBNode[discreHexa[ii]]->getNode();
        };
        mCubicVolume += tetraVolume(vNode[0], vNode[1], vNode[2], vNode[3]);
    };

    return mCubicVolume;
}

// 上位グリッドBNodeへのディレクレ値の分配
//
void CBoundaryHexa::distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();

    double dAveVal;
    uiint iedge, *pnEdgeVert;
    uiint nNumOfEdge=NumberOfEdge::Hexa();
    // 辺のディレクレ値(上位グリッド)
    for(iedge=0; iedge < nNumOfEdge; iedge++){
        pnEdgeVert = pEdgeTree->getHexaLocalNodeNum(iedge);

        dAveVal = 0.0;
        dAveVal += mvBNode[pnEdgeVert[0]]->getValue(dof, mgLevel);
        dAveVal += mvBNode[pnEdgeVert[1]]->getValue(dof, mgLevel);

        dAveVal *= 0.5;//平均値

        if(mnOrder==ElementOrder::Second){
            //2次要素
            mvEdgeBNode[iedge]->setValue(dof, mgLevel,   dAveVal);//カレント・グリッドへDirichlet値をセット

            if(mgLevel!=nMaxMGLevel)
              mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);//上位GridへDirichlet値をセット
        }else{
            //1次要素
            if(mgLevel!=nMaxMGLevel)
              mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);//上位GridへDirichlet値をセット
        }
    };


    CFaceTree *pFaceTree = CFaceTree::Instance();

    uiint iface, *pnFaceVert, ivert;
    uiint nNumOfFace = NumberOfFace::Hexa();

    // 面のディレクレ値(上位グリッド)
    for(iface=0; iface < nNumOfFace; iface++){
        pnFaceVert = pFaceTree->getLocalNodeHexaFace(iface);

        dAveVal=0.0;
        for(ivert=0; ivert < 4; ivert++){
            dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
        };
        dAveVal *= 0.25;//平均値

        if(mgLevel!=nMaxMGLevel)
          mvFaceBNode[iface]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドへディレクレ値をセット
    };


    // 体積中心のディレクレ値(上位グリッド)
    //    &  頂点の値をそのまま上位グリッドへ与える
    double dVertVal;
    dAveVal=0.0;
    uiint nNumOfVert=NumberOfVertex::Hexa();

    for(ivert=0; ivert < nNumOfVert; ivert++){
        dVertVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dAveVal += dVertVal;

        if(mgLevel!=nMaxMGLevel)
          mvBNode[ivert]->setValue(dof, mgLevel+1, dVertVal);//頂点の値を、そのまま上位Gridへセット
    };
    // 体積中心値
    dAveVal *= 0.125;//平均値
    if(mgLevel!=nMaxMGLevel)
        mpVolBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値をセット

////    if(mgLevel > 0){
////        dAveVal *= 0.125;//平均値
////
////        if(mgLevel!=nMaxMGLevel)
////          mpVolBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値をセット
////    }
////    if(mgLevel==0){
////        if(mgLevel!=nMaxMGLevel)
////          mpVolBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
////    }
}

void CBoundaryHexa::replaceEdgeBNode(const uiint& iedge)
{
    uiint nNumOfVert=NumberOfVertex::Hexa();

    CBoundaryNode *pBNode=mvEdgeBNode[iedge];

    // mmBNodeID2Index[pBNode->getID()]= nNumOfVert + iedge;//// <<-- 頂点データ保護のためコメントアウト2011.02.08

    mvBNode[nNumOfVert+iedge]=pBNode;
}


// Refine 後処理 : 辺-面 BNode vectorの解放
// 
void CBoundaryHexa::deleteProgData()
{
    //if(mpElement->getType()==ElementType::Hexa){
    if(mnOrder==ElementOrder::First){
        // Hexa
        vector<CBoundaryNode*>().swap(mvEdgeBNode);// 辺-BNode
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }else{
        // Hexa2
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }
}






