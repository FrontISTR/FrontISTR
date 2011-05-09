//
// BoundaryPrism.cpp
//
//          2010.07.02
//          k.Takeda
#include "BoundaryPrism.h"
#include "ElementProperty.h"
using namespace pmw;


uint CBoundaryPrism::mnElemType = ElementType::Prism;
uint CBoundaryPrism::mNumOfFace = NumberOfFace::Prism();
uint CBoundaryPrism::mNumOfEdge = NumberOfEdge::Prism();
uint CBoundaryPrism::mNumOfNode = NumberOfVertex::Prism();

CBoundaryPrism::CBoundaryPrism()
{
    // Volume形状
    // ----
    // mvbMarkingEdge初期化
    // mvbMarkingFace初期化
    // ----
    uint i;
    // 辺
    mvbMarkingEdge.resize(NumberOfEdge::Prism());
    for(i=0; i < NumberOfEdge::Prism(); i++){
        mvbMarkingEdge[i]=false;
    };
    mvEdgeBNode.resize(NumberOfEdge::Prism());
    mvEdgeNeibVol.resize(NumberOfEdge::Prism());

    // 面
    mvbMarkingFace.resize(NumberOfFace::Prism());
    for(i=0; i < NumberOfFace::Prism(); i++){
        mvbMarkingFace[i]=false;
    };
    mvFaceBNode.resize(NumberOfFace::Prism());
    mvFaceNeibVol.resize(NumberOfFace::Prism());
}
CBoundaryPrism::~CBoundaryPrism()
{
    ;
}


// Volume 形状情報
//
uint CBoundaryPrism::getElemType()
{
    return mnElemType;
}
uint CBoundaryPrism::getNumOfEdge()
{
    return mNumOfEdge;
}

uint CBoundaryPrism::getNumOfFace()
{
    return mNumOfFace;
}

// Edge番号 -> pariBNode
//
PairBNode& CBoundaryPrism::getPairBNode(const uint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uint *pLocalNum;
    
    pLocalNum= pEdgeTree->getPrismLocalNodeNum(iedge);

    uint index1st = pLocalNum[0];
    uint index2nd = pLocalNum[1];

    mPairBNode.first = mvBNode[index1st];
    mPairBNode.second= mvBNode[index2nd];

    return mPairBNode;
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
vector<CBoundaryNode*>& CBoundaryPrism::getFaceCnvNodes(const uint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    CBoundaryNode *pBNode;
    uint  ivert, i, numOfVert;
    uint *pvIndex;

    pvIndex= pFaceTree->getLocalNodePrismFace(iface);
    numOfVert= pFaceTree->getPrismFaceNumOfVert(iface);

    mvFaceCnvNodes.clear();mvFaceCnvNodes.reserve(numOfVert);

    for(i=0; i < numOfVert; i++){
        ivert= pvIndex[i];
        pBNode= mvBNode[ivert];
        mvFaceCnvNodes.push_back(pBNode);
    }

    return mvFaceCnvNodes;
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
        countID++;
        pProgVol->resizeBNode(8);

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
void CBoundaryPrism::distDirichletVal(const uint& dof, const uint& mgLevel)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();

    double dAveVal;
    uint iedge, *pnEdgeVert;

    // 辺のディレクレ値(上位グリッド)
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        pnEdgeVert = pEdgeTree->getPrismLocalNodeNum(iedge);

        dAveVal = 0.0;
        dAveVal += mvBNode[pnEdgeVert[0]]->getValue(dof, mgLevel);
        dAveVal += mvBNode[pnEdgeVert[1]]->getValue(dof, mgLevel);

        dAveVal *= 0.5;//平均値

        mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);//上位GridへDirichlet値をセット
    };

    CFaceTree *pFaceTree = CFaceTree::Instance();

    uint iface, *pnFaceVert, ivert;

    // 面のディレクレ値(上位グリッド)
    for(iface=0; iface < mNumOfFace; iface++){
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

        mvFaceBNode[iface]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドへディレクレ値をセット
    };

    // 体積中心のディレクレ値(上位グリッド)
    //    &  頂点の値をそのまま上位グリッドへ与える
    double dVertVal;
    dAveVal=0.0;
    for(ivert=0; ivert < mNumOfNode; ivert++){
        dVertVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dAveVal += dVertVal;

        mvBNode[ivert]->setValue(dof, mgLevel+1, dVertVal);//頂点の値は、そのまま上位Gridへセット
    };
    if(mgLevel > 0){
        dAveVal /= 6.0;//平均値
        mpVolBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値をセット
    }
    if(mgLevel==0){
        mpVolBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
    }
}





