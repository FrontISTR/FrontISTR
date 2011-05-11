//
// BoundaryTetra.cpp
//
//              2010.07.02
//              k.Takeda
#include "BoundaryTetra.h"
#include "ElementProperty.h"
using namespace pmw;



//uint CBoundaryTetra::mnElemType = ElementType::Tetra;
//uint CBoundaryTetra::mNumOfFace = NumberOfFace::Tetra();
//uint CBoundaryTetra::mNumOfEdge = NumberOfEdge::Tetra();
//uint CBoundaryTetra::mNumOfNode = NumberOfVertex::Tetra();


CBoundaryTetra::CBoundaryTetra()
{
    // Volume形状
    // ----
    // mvbMarkingEdge初期化
    // mvbMarkingFace初期化
    // ----
    uiint i;
    // 辺
    mvbMarkingEdge = new bool[NumberOfEdge::Tetra()];
    for(i=0; i < NumberOfEdge::Tetra(); i++){
        mvbMarkingEdge[i]=false;
    };
    mvEdgeBNode.resize(NumberOfEdge::Tetra());
    mvEdgeNeibVol.resize(NumberOfEdge::Tetra());

    // 面
    mvbMarkingFace = new bool[NumberOfFace::Tetra()];
    for(i=0; i < NumberOfFace::Tetra(); i++){
        mvbMarkingFace[i]=false;
    };
    mvFaceBNode.resize(NumberOfFace::Tetra());
    mvFaceNeibVol.resize(NumberOfFace::Tetra());

    for(i=0; i < NumberOfFace::Tetra(); i++){
        mvFaceBNode[i]= NULL;
    }
    mpVolBNode= NULL;
}

CBoundaryTetra::~CBoundaryTetra()
{
    ;
}



// Volume 形状情報
//
uiint CBoundaryTetra::getElemType()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();

    switch(mnOrder){
        case(ElementOrder::First):
            return ElementType::Tetra;
        case(ElementOrder::Second):
            return ElementType::Tetra2;
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryTetra::getElemType, invalid mnOrder");
            return ElementType::Limit;
    }
}
uiint CBoundaryTetra::getNumOfEdge()
{
    return NumberOfEdge::Tetra();
}

uiint CBoundaryTetra::getNumOfFace()
{
    return NumberOfFace::Tetra();
}

uiint CBoundaryTetra::getNumOfNode()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();

    switch(mnOrder){
        case(ElementOrder::First):
            return NumberOfNode::Tetra();
        case(ElementOrder::Second):
            return NumberOfNode::Tetra2();
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryTetra::getNumOfNode, invalid mnOrder");
            return NumberOfNode::Default();
    }
}

uiint CBoundaryTetra::getNumOfVert()
{
    return NumberOfVertex::Tetra();
}

void CBoundaryTetra::setOrder(const uiint& order)
{
    mnOrder = order;

    switch(mnOrder){
    case(ElementOrder::First):
        mvBNode.resize(NumberOfNode::Tetra());
        break;

    case(ElementOrder::Second):
        mvBNode.resize(NumberOfNode::Tetra2());
        break;
    }
}

// Edge番号 -> pariBNode
//
PairBNode CBoundaryTetra::getPairBNode(const uiint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uiint *pLocalNum;
    
    pLocalNum= pEdgeTree->getTetraLocalNodeNum(iedge);
     
    uiint index1st = pLocalNum[0];
    uiint index2nd = pLocalNum[1];

    PairBNode pairBNode;
    pairBNode.first = mvBNode[index1st];
    pairBNode.second= mvBNode[index2nd];

    return pairBNode;
}


// pairBNode -> Edge番号
//
uiint& CBoundaryTetra::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uiint ivert= getVertIndex(pairBNode.first);
    uiint jvert= getVertIndex(pairBNode.second);
    
    return pEdgeTree->getTetraEdgeIndex(ivert, jvert);
}


// 面を構成するBNode配列
//
vector<CBoundaryNode*> CBoundaryTetra::getFaceCnvNodes(const uiint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    CBoundaryNode *pBNode;
    uiint  ivert, i, numOfVert;
    uiint *pvIndex;
    
    pvIndex= pFaceTree->getLocalNodeTetraFace(iface);
    numOfVert= pFaceTree->getTetraFaceNumOfVert(iface);

    vector<CBoundaryNode*> vFaceCnvNodes;
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
uiint& CBoundaryTetra::getFaceID(vector<CBoundaryNode*>& vBNode)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();

    vuint vLocalVert;
    CBoundaryNode *pBNode;
    uiint ibnode, numOfBNode= vBNode.size();

    vLocalVert.reserve(3);
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= vBNode[ibnode];
        vLocalVert.push_back(getVertIndex(pBNode));
    };

    return pFaceTree->getTetraFaceIndex2(vLocalVert);
}

// 辺を構成する頂点番号
uiint* CBoundaryTetra::getLocalNode_Edge(const uiint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    return pEdgeTree->getTetraLocalNodeNum(iedge);
}
// 面を構成する頂点番号
uiint* CBoundaryTetra::getLocalNode_Face(const uiint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();

    return pFaceTree->getLocalNodeTetraFace(iface);
}




// Refine Volume再分割
// ----
void CBoundaryTetra::refine(uiint& countID, const vuint& vDOF)
{
    CBoundaryVolume *pProgVol;
    CElement   *pProgElem;
    CNode      *pNode;
    CBoundaryNode *pBNode;
    uiint iprog, numOfProg;
    uiint nProgPos;

    double progCubicVol;//体積
    double coef;

    numOfProg= 4;
    mvProgVolume.reserve(numOfProg);

    for(iprog=0; iprog < numOfProg; iprog++){

        pProgVol = new CBoundaryHexa;// <<<<<<<<<<<<<<<<<<<<<<< new

        pProgVol->setID(countID);
        pProgVol->setOrder(mnOrder);//resizeBNodeも含む
        countID++;

        //pProgVol->resizeBNode(8);

        mvProgVolume.push_back(pProgVol);

        nProgPos = dividTetra(iprog, pProgVol);//子供がぶら下がっている場所::nProgPosはBNodeの頂点番号(配列番号)

        //子要素を取得
        pBNode= mvBNode[nProgPos];//
        pNode = pBNode->getNode();
        pProgElem= mpElement->getProgElem_NodeID(pNode->getID());// <<<<<<<<<<<<<<<<<<<<< Nodeにぶらさがっている
        // 子Volに子要素をセット
        pProgVol->setElement(pProgElem);
        pProgVol->setElementID(pProgElem->getID());

        // 体積比による境界値の分配
        // ----
        progCubicVol= pProgVol->calcVolume();
        coef = progCubicVol/mCubicVolume;

        distValue(pProgVol, coef, vDOF);

    };// iprog loop_end

}

// BoundaryVolumeの体積
// ----
double& CBoundaryTetra::calcVolume()
{
    uiint  i;
    CNode* vNode[4];

    for(i=0; i< 4; i++){
        vNode[i]= mvBNode[i]->getNode();
    };
    mCubicVolume = tetraVolume(vNode[0], vNode[1], vNode[2], vNode[3]);

    return mCubicVolume;
}


// 上位グリッドBNodeへのディレクレ値の分配
//
void CBoundaryTetra::distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();

    double dAveVal;
    uiint iedge, *pnEdgeVert;
    uiint nNumOfEdge=NumberOfEdge::Tetra();

    // 辺のディレクレ値(上位グリッド)
    for(iedge=0; iedge < nNumOfEdge; iedge++){
        pnEdgeVert = pEdgeTree->getTetraLocalNodeNum(iedge);

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
    uiint nNumOfFace=NumberOfFace::Tetra();

    // 面のディレクレ値(上位グリッド)
    for(iface=0; iface < nNumOfFace; iface++){
        pnFaceVert = pFaceTree->getLocalNodeTetraFace(iface);

        dAveVal=0.0;
        for(ivert=0; ivert < 3; ivert++){
            dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
        };
        dAveVal /= 3.0;//平均値

        if(mgLevel!=nMaxMGLevel)
          mvFaceBNode[iface]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドへディレクレ値をセット
        
    };

    // 体積中心のディレクレ値(上位グリッド)
    //    &  頂点の値をそのまま上位グリッドへ与える
    double dVertVal;
    dAveVal=0.0;
    uiint nNumOfVert=NumberOfVertex::Tetra();

    for(ivert=0; ivert < nNumOfVert; ivert++){
        dVertVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dAveVal += dVertVal;

        if(mgLevel!=nMaxMGLevel)
          mvBNode[ivert]->setValue(dof, mgLevel+1, dVertVal);//頂点の値を、そのまま上位Gridへセット
    };
    // 体積中心値
    dAveVal *= 0.25;//平均値
    if(mgLevel!=nMaxMGLevel)
      mpVolBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値をセット

////    if(mgLevel > 0){
////        dAveVal *= 0.25;//平均値
////
////        if(mgLevel!=nMaxMGLevel)
////          mpVolBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値をセット
////    }
////    if(mgLevel==0){
////        if(mgLevel!=nMaxMGLevel)
////          mpVolBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
////    }
}

void CBoundaryTetra::replaceEdgeBNode(const uiint& iedge)
{
    uiint nNumOfVert=NumberOfVertex::Tetra();

    CBoundaryNode *pBNode=mvEdgeBNode[iedge];

    // mmBNodeID2Index[pBNode->getID()] = nNumOfVert + iedge;//// <<-- 頂点データ保護のためコメントアウト2011.02.08

    mvBNode[nNumOfVert+iedge]=pBNode;
}

// Refine 後処理 : 辺-面 BNode vectorの解放
//
void CBoundaryTetra::deleteProgData()
{
    //if(mpElement->getType()==ElementType::Tetra){
    if(mnOrder==ElementOrder::First){
        // Tetra
        vector<CBoundaryNode*>().swap(mvEdgeBNode);// 辺-BNode
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }else{
        // Tetra2
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }
}




