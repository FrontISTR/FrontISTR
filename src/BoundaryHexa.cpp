//
// BoundaryHexa.cpp
//
//          2010.07.02
//          k.Takeda
#include "BoundaryHexa.h"
#include "ElementProperty.h"
using namespace pmw;


uint CBoundaryHexa::mnElemType = ElementType::Hexa;
uint CBoundaryHexa::mNumOfFace = NumberOfFace::Hexa();
uint CBoundaryHexa::mNumOfEdge = NumberOfEdge::Hexa();
uint CBoundaryHexa::mNumOfNode = NumberOfVertex::Hexa();

CBoundaryHexa::CBoundaryHexa()
{
    // Volume形状
    // ----
    // mvbMarkingEdge初期化
    // mvbMarkingFace初期化
    // ----
    uint i;
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
}
CBoundaryHexa::~CBoundaryHexa()
{
    ;
}


// Volume 形状情報
//
uint CBoundaryHexa::getElemType()
{
    return mnElemType;
}
uint CBoundaryHexa::getNumOfEdge()
{
    return mNumOfEdge;
}
uint CBoundaryHexa::getNumOfFace()
{
    return mNumOfFace;
}

// Edge番号 -> pariBNode
//
PairBNode CBoundaryHexa::getPairBNode(const uint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uint *pLocalNum;
    
    pLocalNum= pEdgeTree->getHexaLocalNodeNum(iedge);

    uint index1st = pLocalNum[0];
    uint index2nd = pLocalNum[1];

    PairBNode pairBNode;
    pairBNode.first = mvBNode[index1st];
    pairBNode.second= mvBNode[index2nd];

    return pairBNode;
}


// pairBNode -> Edge番号
//
uint& CBoundaryHexa::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uint ivert= getVertIndex(pairBNode.first);
    uint jvert= getVertIndex(pairBNode.second);

    return pEdgeTree->getHexaEdgeIndex(ivert,jvert);
}


// 面を構成するBNode配列
//
vector<CBoundaryNode*> CBoundaryHexa::getFaceCnvNodes(const uint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    CBoundaryNode *pBNode;
    uint  ivert, i, numOfVert;
    uint *pvIndex;
    
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
uint& CBoundaryHexa::getFaceID(vector<CBoundaryNode*>& vBNode)
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

    return pFaceTree->getHexaFaceIndex2(vLocalVert);
}

// 辺を構成する頂点番号
uint* CBoundaryHexa::getLocalNode_Edge(const uint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    return pEdgeTree->getHexaLocalNodeNum(iedge);
}
// 面を構成する頂点番号
uint* CBoundaryHexa::getLocalNode_Face(const uint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();

    return pFaceTree->getLocalNodeHexaFace(iface);
}



// Refine Volume再分割
// ----
void CBoundaryHexa::refine(uint& countID, const vuint& vDOF)
{
    CBoundaryVolume *pProgVol;
    CElement   *pProgElem;
    CNode      *pNode;
    CBoundaryNode *pBNode;
    uint iprog, numOfProg;
    uint nProgPos;
    double progCubicVol;//体積
    double coef;

    numOfProg= 8;
    mvProgVolume.reserve(numOfProg);

    for(iprog=0; iprog < numOfProg; iprog++){

        pProgVol = new CBoundaryHexa;// <<<<<<<<<<<<<<<< new

        pProgVol->setID(countID);
        countID++;
        pProgVol->resizeBNode(8);

        mvProgVolume.push_back(pProgVol);

        nProgPos = dividHexa(iprog, pProgVol);//子供がぶら下がっている場所::nProgPosはBNodeの頂点番号(配列番号)

        //子要素を取得
        pBNode= mvBNode[nProgPos];
        pNode = pBNode->getNode();
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
    uint* discreHexa;
    uint  i,ii;
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
void CBoundaryHexa::distDirichletVal(const uint& dof, const uint& mgLevel)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();
    
    double dAveVal;
    uint iedge, *pnEdgeVert;
    
    // 辺のディレクレ値(上位グリッド)
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        pnEdgeVert = pEdgeTree->getHexaLocalNodeNum(iedge);
        
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
        pnFaceVert = pFaceTree->getLocalNodeHexaFace(iface);
        
        dAveVal=0.0;
        for(ivert=0; ivert < 4; ivert++){
            dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
        };
        dAveVal *= 0.25;//平均値
        
        mvFaceBNode[iface]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドへディレクレ値をセット
    };

    
    // 体積中心のディレクレ値(上位グリッド)
    //    &  頂点の値をそのまま上位グリッドへ与える
    double dVertVal;
    dAveVal=0.0;
    for(ivert=0; ivert < mNumOfNode; ivert++){
        dVertVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dAveVal += dVertVal;

        mvBNode[ivert]->setValue(dof, mgLevel+1, dVertVal);//頂点の値を、そのまま上位Gridへセット
    };
    if(mgLevel > 0){
        dAveVal *= 0.125;//平均値
        mpVolBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値をセット
    }
    if(mgLevel==0){
        mpVolBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
    }
}


// Refine 後処理 : 辺-面 BNode vectorの解放
// 
void CBoundaryHexa::deleteProgData()
{
    if(mpElement->getType()==ElementType::Hexa){
        vector<CBoundaryNode*>().swap(mvEdgeBNode);// 辺-BNode
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }else{
        // Hexa2
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }
}






