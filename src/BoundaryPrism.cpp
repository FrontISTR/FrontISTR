//
// BoundaryPrism.cpp
//
//          2010.07.02
//          k.Takeda
#include "BoundaryPrism.h"
#include "ElementProperty.h"
using namespace pmw;


//uint CBoundaryPrism::mnElemType = ElementType::Prism;
//uint CBoundaryPrism::mNumOfFace = NumberOfFace::Prism();
//uint CBoundaryPrism::mNumOfEdge = NumberOfEdge::Prism();
//uint CBoundaryPrism::mNumOfNode = NumberOfVertex::Prism();

CBoundaryPrism::CBoundaryPrism()
{
    // Volume形状
    // ----
    // mvbMarkingEdge初期化
    // mvbMarkingFace初期化
    // ----
    uint i;
    // 辺
    mvbMarkingEdge = new bool[NumberOfEdge::Prism()];
    for(i=0; i < NumberOfEdge::Prism(); i++){
        mvbMarkingEdge[i]=false;
    };
    mvEdgeBNode.resize(NumberOfEdge::Prism());
    mvEdgeNeibVol.resize(NumberOfEdge::Prism());

    // 面
    mvbMarkingFace = new bool[NumberOfFace::Prism()];
    for(i=0; i < NumberOfFace::Prism(); i++){
        mvbMarkingFace[i]=false;
    };
    mvFaceBNode.resize(NumberOfFace::Prism());
    mvFaceNeibVol.resize(NumberOfFace::Prism());

    for(i=0; i < NumberOfFace::Prism(); i++){
        mvFaceBNode[i]= NULL;
    }
    mpVolBNode= NULL;
}
CBoundaryPrism::~CBoundaryPrism()
{
    ;
}


// Volume 形状情報
//
uint CBoundaryPrism::getElemType()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();

    switch(mnOrder){
        case(ElementOrder::First):
            return ElementType::Prism;
        case(ElementOrder::Second):
            return ElementType::Prism2;
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryPrism::getElemType, invalid mnOrder");
            return ElementType::Limit;
    }
}
uint CBoundaryPrism::getNumOfEdge()
{
    return NumberOfEdge::Prism();
}

uint CBoundaryPrism::getNumOfFace()
{
    return NumberOfFace::Prism();
}

uint CBoundaryPrism::getNumOfNode()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();

    switch(mnOrder){
        case(ElementOrder::First):
            return NumberOfNode::Prism();
        case(ElementOrder::Second):
            return NumberOfNode::Prism2();
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryPrism::getNumOfNode, invalid mnOrder");
            return NumberOfNode::Default();
    }
}

uint CBoundaryPrism::getNumOfVert()
{
    return NumberOfVertex::Prism();;
}

void CBoundaryPrism::setOrder(const uint& order)
{
    mnOrder = order;

    switch(mnOrder){
    case(ElementOrder::First):
        mvBNode.resize(NumberOfNode::Prism());
        break;

    case(ElementOrder::Second):
        mvBNode.resize(NumberOfNode::Prism2());
        break;
    }
}

// Edge番号 -> pariBNode
//
PairBNode CBoundaryPrism::getPairBNode(const uint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uint *pLocalNum;
    
    pLocalNum= pEdgeTree->getPrismLocalNodeNum(iedge);

    uint index1st = pLocalNum[0];
    uint index2nd = pLocalNum[1];

    PairBNode pairBNode;
    pairBNode.first = mvBNode[index1st];
    pairBNode.second= mvBNode[index2nd];

    return pairBNode;
}


// pairBNode -> Edge番号
//
uint& CBoundaryPrism::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uint ivert= getVertIndex(pairBNode.first);
    uint jvert= getVertIndex(pairBNode.second);

    return pEdgeTree->getPrismEdgeIndex(ivert, jvert);
}


// 面を構成するBNode配列
//
vector<CBoundaryNode*> CBoundaryPrism::getFaceCnvNodes(const uint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    CBoundaryNode *pBNode;
    uint  ivert, i, numOfVert;
    uint *pvIndex;

    pvIndex= pFaceTree->getLocalNodePrismFace(iface);
    numOfVert= pFaceTree->getPrismFaceNumOfVert(iface);

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
uint& CBoundaryPrism::getFaceID(vector<CBoundaryNode*>& vBNode)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();

    vuint vLocalVert;
    CBoundaryNode *pBNode;
    uint ibnode, numOfBNode= vBNode.size();

    vLocalVert.reserve(4);
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= vBNode[ibnode];
        vLocalVert.push_back(getVertIndex(pBNode));
    };

    return pFaceTree->getPrismFaceIndex2(vLocalVert);
}


// 辺を構成する頂点番号
uint* CBoundaryPrism::getLocalNode_Edge(const uint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    return pEdgeTree->getPrismLocalNodeNum(iedge);
}
// 面を構成する頂点番号
uint* CBoundaryPrism::getLocalNode_Face(const uint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();

    return pFaceTree->getLocalNodePrismFace(iface);
}



// Refine Volume再分割
// ----
void CBoundaryPrism::refine(uint& countID, const vuint& vDOF)
{
    CBoundaryVolume *pProgVol;
    CElement   *pProgElem;
    CNode      *pNode;
    CBoundaryNode *pBNode;
    uint iprog, numOfProg;
    uint nProgPos;

    double progCubicVol;//体積
    double coef;

    numOfProg= 6;
    mvProgVolume.reserve(numOfProg);

    for(iprog=0; iprog < numOfProg; iprog++){

        pProgVol = new CBoundaryHexa;// <<<<<<<<<<<<<<<<< new

        pProgVol->setID(countID);
        pProgVol->setOrder(mnOrder);//resizeBNodeも含む
        countID++;
        //pProgVol->resizeBNode(8);

        mvProgVolume.push_back(pProgVol);

        nProgPos = dividPrism(iprog, pProgVol);//子供がぶら下がっている場所::nProgPosはBNodeの頂点番号(配列番号)

        //子要素を取得
        pBNode= mvBNode[nProgPos];
        pNode = pBNode->getNode();
        pProgElem= mpElement->getProgElem_NodeID(pNode->getID());// <<<<<<<<<<<<< Nodeにぶらさがっている
        // 子Volに子要素をセット
        pProgVol->setElement(pProgElem);
        pProgVol->setElementID(pProgElem->getID());


        // 体積比による境界値の分配
        progCubicVol= pProgVol->calcVolume();
        coef = progCubicVol/mCubicVolume;

        distValue(pProgVol, coef, vDOF);

    };// iprog loop_end

}

// BoundaryVolumeの体積
// ----
double& CBoundaryPrism::calcVolume()
{
    CDiscreteVolume *pDiscrete= CDiscreteVolume::Instance();
    uint* discrePrism;
    uint  i,ii;
    CNode* vNode[4];

    mCubicVolume= 0.0;

    //3個のTetraの体積を加算
    for(i=0; i< 3; i++){
        discrePrism= pDiscrete->getPrismDiscrete(i);
        for(ii=0; ii< 4; ii++){
            vNode[ii]= mvBNode[discrePrism[ii]]->getNode();
        };
        mCubicVolume += tetraVolume(vNode[0], vNode[1], vNode[2], vNode[3]);
    };

    return mCubicVolume;
}


// 上位グリッドBNodeへのディレクレ値の分配
//
void CBoundaryPrism::distDirichletVal(const uint& dof, const uint& mgLevel, const uint& nMaxMGLevel)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();

    double dAveVal;
    uint iedge, *pnEdgeVert;
    uint nNumOfEdge=NumberOfEdge::Prism();

    // 辺のディレクレ値(上位グリッド)
    for(iedge=0; iedge < nNumOfEdge; iedge++){
        pnEdgeVert = pEdgeTree->getPrismLocalNodeNum(iedge);

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

    uint iface, *pnFaceVert, ivert;
    uint nNumOfFace=NumberOfFace::Prism();

    // 面のディレクレ値(上位グリッド)
    for(iface=0; iface < nNumOfFace; iface++){
        pnFaceVert = pFaceTree->getLocalNodePrismFace(iface);

        dAveVal=0.0;
        if(iface > 1){
            for(ivert=0; ivert < 4; ivert++){
                dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
            };
            dAveVal *= 0.25;//平均値
        }else{
            for(ivert=0; ivert < 3; ivert++){
                dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
            };
            dAveVal /= 3.0;//平均値
        }

        if(mgLevel!=nMaxMGLevel)
          mvFaceBNode[iface]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドへディレクレ値をセット
    };

    // 体積中心のディレクレ値(上位グリッド)
    //    &  頂点の値をそのまま上位グリッドへ与える
    double dVertVal;
    dAveVal=0.0;
    uint nNumOfVert=NumberOfVertex::Prism();

    for(ivert=0; ivert < nNumOfVert; ivert++){
        dVertVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dAveVal += dVertVal;

        if(mgLevel!=nMaxMGLevel)
          mvBNode[ivert]->setValue(dof, mgLevel+1, dVertVal);//頂点の値は、そのまま上位Gridへセット
    };
    if(mgLevel > 0){
        dAveVal /= 6.0;//平均値
        if(mgLevel!=nMaxMGLevel)
          mpVolBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値をセット
    }
    if(mgLevel==0){
        if(mgLevel!=nMaxMGLevel)
          mpVolBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
    }
}

void CBoundaryPrism::replaceEdgeBNode(const uint& iedge)
{
    uint nNumOfVert=NumberOfVertex::Prism();

    CBoundaryNode *pBNode=mvEdgeBNode[iedge];

    // mmBNodeID2Index[pBNode->getID()] = nNumOfVert + iedge;//// <<-- 頂点データ保護のためコメントアウト2011.02.08

    mvBNode[nNumOfVert+iedge]=pBNode;
}

// Refine 後処理 : 辺-面 BNode vectorの解放
//
void CBoundaryPrism::deleteProgData()
{
    //if(mpElement->getType()==ElementType::Prism){
    if(mnOrder==ElementOrder::First){
        // Prism
        vector<CBoundaryNode*>().swap(mvEdgeBNode);// 辺-BNode
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }else{
        // Prism2
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }
}



