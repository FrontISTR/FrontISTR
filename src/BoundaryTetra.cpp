//
// BoundaryTetra.cpp
//
//              2010.07.02
//              k.Takeda
#include "BoundaryTetra.h"
#include "ElementProperty.h"
using namespace pmw;



uint CBoundaryTetra::mnElemType = ElementType::Tetra;
uint CBoundaryTetra::mNumOfFace = NumberOfFace::Tetra();
uint CBoundaryTetra::mNumOfEdge = NumberOfEdge::Tetra();
uint CBoundaryTetra::mNumOfNode = NumberOfVertex::Tetra();


CBoundaryTetra::CBoundaryTetra()
{
    // Volume形状
    // ----
    // mvbMarkingEdge初期化
    // mvbMarkingFace初期化
    // ----
    uint i;
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
}

CBoundaryTetra::~CBoundaryTetra()
{
    ;
}



// Volume 形状情報
//
uint CBoundaryTetra::getElemType()
{
    return mnElemType;
}
uint CBoundaryTetra::getNumOfEdge()
{
    return mNumOfEdge;
}

uint CBoundaryTetra::getNumOfFace()
{
    return mNumOfFace;
}

// Edge番号 -> pariBNode
//
PairBNode CBoundaryTetra::getPairBNode(const uint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uint *pLocalNum;
    
    pLocalNum= pEdgeTree->getTetraLocalNodeNum(iedge);
     
    uint index1st = pLocalNum[0];
    uint index2nd = pLocalNum[1];

    PairBNode pairBNode;
    pairBNode.first = mvBNode[index1st];
    pairBNode.second= mvBNode[index2nd];

    return pairBNode;
}


// pairBNode -> Edge番号
//
uint& CBoundaryTetra::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uint ivert= getVertIndex(pairBNode.first);
    uint jvert= getVertIndex(pairBNode.second);
    
    return pEdgeTree->getTetraEdgeIndex(ivert, jvert);
}


// 面を構成するBNode配列
//
vector<CBoundaryNode*> CBoundaryTetra::getFaceCnvNodes(const uint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    CBoundaryNode *pBNode;
    uint  ivert, i, numOfVert;
    uint *pvIndex;
    
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
uint& CBoundaryTetra::getFaceID(vector<CBoundaryNode*>& vBNode)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();

    vuint vLocalVert;
    CBoundaryNode *pBNode;
    uint ibnode, numOfBNode= vBNode.size();

    vLocalVert.reserve(3);
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= vBNode[ibnode];
        vLocalVert.push_back(getVertIndex(pBNode));
    };

    return pFaceTree->getTetraFaceIndex2(vLocalVert);
}

// 辺を構成する頂点番号
uint* CBoundaryTetra::getLocalNode_Edge(const uint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    return pEdgeTree->getTetraLocalNodeNum(iedge);
}
// 面を構成する頂点番号
uint* CBoundaryTetra::getLocalNode_Face(const uint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();

    return pFaceTree->getLocalNodeTetraFace(iface);
}




// Refine Volume再分割
// ----
void CBoundaryTetra::refine(uint& countID, const vuint& vDOF)
{
    CBoundaryVolume *pProgVol;
    CElement   *pProgElem;
    CNode      *pNode;
    CBoundaryNode *pBNode;
    uint iprog, numOfProg;
    uint nProgPos;

    double progCubicVol;//体積
    double coef;

    numOfProg= 4;
    mvProgVolume.reserve(numOfProg);

    for(iprog=0; iprog < numOfProg; iprog++){

        pProgVol = new CBoundaryHexa;// <<<<<<<<<<<<<<<<<<<<<<< new

        pProgVol->setID(countID);
        countID++;

        pProgVol->resizeBNode(8);

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
    uint  i;
    CNode* vNode[4];

    for(i=0; i< 4; i++){
        vNode[i]= mvBNode[i]->getNode();
    };
    mCubicVolume = tetraVolume(vNode[0], vNode[1], vNode[2], vNode[3]);

    return mCubicVolume;
}


// 上位グリッドBNodeへのディレクレ値の分配
//
void CBoundaryTetra::distDirichletVal(const uint& dof, const uint& mgLevel)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();

    double dAveVal;
    uint iedge, *pnEdgeVert;

    // 辺のディレクレ値(上位グリッド)
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        pnEdgeVert = pEdgeTree->getTetraLocalNodeNum(iedge);

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
        pnFaceVert = pFaceTree->getLocalNodeTetraFace(iface);

        dAveVal=0.0;
        for(ivert=0; ivert < 3; ivert++){
            dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
        };
        dAveVal /= 3.0;//平均値

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
        dAveVal *= 0.25;//平均値
        mpVolBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値をセット
    }
    if(mgLevel==0){
        mpVolBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
    }
}


// Refine 後処理 : 辺-面 BNode vectorの解放
//
void CBoundaryTetra::deleteProgData()
{
    if(mpElement->getType()==ElementType::Tetra){
        vector<CBoundaryNode*>().swap(mvEdgeBNode);// 辺-BNode
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }else{
        // Tetra2
        vector<CBoundaryNode*>().swap(mvFaceBNode);// 面-BNode
    }
}




