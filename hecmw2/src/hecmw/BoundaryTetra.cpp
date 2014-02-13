/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryTetra.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "BoundaryTetra.h"
#include "ElementProperty.h"
using namespace pmw;
CBoundaryTetra::CBoundaryTetra()
{
    uiint i;
    mvbMarkingEdge = new bool[NumberOfEdge::Tetra()];
    for(i=0; i < NumberOfEdge::Tetra(); i++) {
        mvbMarkingEdge[i]=false;
    };
    mvEdgeBNode.resize(NumberOfEdge::Tetra());
    mvEdgeNeibVol.resize(NumberOfEdge::Tetra());
    mvbMarkingFace = new bool[NumberOfFace::Tetra()];
    for(i=0; i < NumberOfFace::Tetra(); i++) {
        mvbMarkingFace[i]=false;
    };
    mvFaceBNode.resize(NumberOfFace::Tetra());
    mvFaceNeibVol.resize(NumberOfFace::Tetra());
    for(i=0; i < NumberOfFace::Tetra(); i++) {
        mvFaceBNode[i]= NULL;
    }
    mpVolBNode= NULL;
}
CBoundaryTetra::~CBoundaryTetra()
{
    ;
}
uiint CBoundaryTetra::getElemType()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    switch(mnOrder) {
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
    switch(mnOrder) {
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
    switch(mnOrder) {
    case(ElementOrder::First):
        mvBNode.resize(NumberOfNode::Tetra());
        break;
    case(ElementOrder::Second):
        mvBNode.resize(NumberOfNode::Tetra2());
        break;
    }
}
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
uiint& CBoundaryTetra::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uiint ivert= getVertIndex(pairBNode.first);
    uiint jvert= getVertIndex(pairBNode.second);
    return pEdgeTree->getTetraEdgeIndex(ivert, jvert);
}
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
    for(i=0; i < numOfVert; i++) {
        ivert= pvIndex[i];
        pBNode= mvBNode[ivert];
        vFaceCnvNodes.push_back(pBNode);
    }
    return vFaceCnvNodes;
}
uiint& CBoundaryTetra::getFaceID(vector<CBoundaryNode*>& vBNode)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    vuint vLocalVert;
    CBoundaryNode *pBNode;
    uiint ibnode, numOfBNode= vBNode.size();
    vLocalVert.reserve(3);
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= vBNode[ibnode];
        vLocalVert.push_back(getVertIndex(pBNode));
    };
    return pFaceTree->getTetraFaceIndex2(vLocalVert);
}
uiint* CBoundaryTetra::getLocalNode_Edge(const uiint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getTetraLocalNodeNum(iedge);
}
uiint* CBoundaryTetra::getLocalNode_Face(const uiint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    return pFaceTree->getLocalNodeTetraFace(iface);
}
void CBoundaryTetra::refine(uiint& countID, const vuint& vDOF)
{
    CBoundaryVolume *pProgVol;
    CElement   *pProgElem;
    CNode      *pNode;
    CBoundaryNode *pBNode;
    uiint iprog, numOfProg;
    uiint nProgPos;
    double progCubicVol;
    double coef;
    numOfProg= 4;
    mvProgVolume.reserve(numOfProg);
    for(iprog=0; iprog < numOfProg; iprog++) {
        pProgVol = new CBoundaryHexa;
        pProgVol->setID(countID);
        pProgVol->setOrder(mnOrder);
        countID++;
        mvProgVolume.push_back(pProgVol);
        nProgPos = dividTetra(iprog, pProgVol);
        pBNode= mvBNode[nProgPos];
        pNode = pBNode->getNode();
        pProgElem= mpElement->getProgElem_NodeID(pNode->getID());
        pProgVol->setElement(pProgElem);
        pProgVol->setElementID(pProgElem->getID());
        progCubicVol= pProgVol->calcVolume();
        coef = progCubicVol/mCubicVolume;
        distValue(pProgVol, coef, vDOF);
    };
}
double& CBoundaryTetra::calcVolume()
{
    uiint  i;
    CNode* vNode[4];
    for(i=0; i< 4; i++) {
        vNode[i]= mvBNode[i]->getNode();
    };
    mCubicVolume = tetraVolume(vNode[0], vNode[1], vNode[2], vNode[3]);
    return mCubicVolume;
}
void CBoundaryTetra::distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();
    double dAveVal, calcVal;
    uiint iedge, *pnEdgeVert;
    uiint nNumOfEdge=NumberOfEdge::Tetra();
    // 辺中心
    for(iedge=0; iedge < nNumOfEdge; iedge++) {
        pnEdgeVert = pEdgeTree->getTetraLocalNodeNum(iedge);
        dAveVal = 0.0;
        //dAveVal += mvBNode[pnEdgeVert[0]]->getValue(dof, mgLevel);
        //dAveVal += mvBNode[pnEdgeVert[1]]->getValue(dof, mgLevel);
        dAveVal += mvBNode[pnEdgeVert[0]]->getEntValue(dof, mgLevel);
        dAveVal += mvBNode[pnEdgeVert[1]]->getEntValue(dof, mgLevel);
        dAveVal *= 0.5;
        if(mnOrder==ElementOrder::Second) {
            //mvEdgeBNode[iedge]->setValue(dof, mgLevel,   dAveVal);
            mvEdgeBNode[iedge]->setEntValue(dof, mgLevel,   dAveVal);

            calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
            mvEdgeBNode[iedge]->setValue(dof, mgLevel,  calcVal);//-----Level
        }
        if(mgLevel!=nMaxMGLevel) {
            //mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);
            mvEdgeBNode[iedge]->setEntValue(dof, mgLevel+1, dAveVal);

            calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
            mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, calcVal);//---------Level+1
        }

    };
    // 面中心
    CFaceTree *pFaceTree = CFaceTree::Instance();
    uiint iface, *pnFaceVert, ivert;
    uiint nNumOfFace=NumberOfFace::Tetra();
    for(iface=0; iface < nNumOfFace; iface++) {
        pnFaceVert = pFaceTree->getLocalNodeTetraFace(iface);
        dAveVal=0.0;
        for(ivert=0; ivert < 3; ivert++) {
            //dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
            dAveVal += mvBNode[pnFaceVert[ivert]]->getEntValue(dof, mgLevel);
        };
        dAveVal /= 3.0;
        if(mgLevel!=nMaxMGLevel) {
            //mvFaceBNode[iface]->setValue(dof, mgLevel+1, dAveVal);
            mvFaceBNode[iface]->setEntValue(dof, mgLevel+1, dAveVal);

            calcVal= getCalcValue( mvFaceBNode[iface], dAveVal, dof, vPolandDOF, mPoland );
            mvFaceBNode[iface]->setValue(dof, mgLevel+1, calcVal);//----------Level+1
        }
    };
    // 頂点
    double dVertVal;
    dAveVal=0.0;
    uiint nNumOfVert=NumberOfVertex::Tetra();
    for(ivert=0; ivert < nNumOfVert; ivert++) {
        //dVertVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dVertVal = mvBNode[ivert]->getEntValue(dof, mgLevel);
        dAveVal += dVertVal;
        if(mgLevel!=nMaxMGLevel) {
            //mvBNode[ivert]->setValue(dof, mgLevel+1, dVertVal);
            mvBNode[ivert]->setEntValue(dof, mgLevel+1, dVertVal);

            calcVal= getCalcValue( mvBNode[ivert], dVertVal, dof, vPolandDOF, mPoland );
            mvBNode[ivert]->setValue(dof, mgLevel+1, calcVal);//-------------Level+1
        }
    };
    // 要素中心
    dAveVal *= 0.25;
    if(mgLevel!=nMaxMGLevel) {
        //mpVolBNode->setValue(dof, mgLevel+1, dAveVal);
        mpVolBNode->setEntValue(dof, mgLevel+1, dAveVal);

        calcVal= getCalcValue( mpVolBNode, dAveVal, dof, vPolandDOF, mPoland );
        mpVolBNode->setValue(dof, mgLevel+1, calcVal);//-------------Level+1
    }
}
void CBoundaryTetra::replaceEdgeBNode(const uiint& iedge)
{
    uiint nNumOfVert=NumberOfVertex::Tetra();
    CBoundaryNode *pBNode=mvEdgeBNode[iedge];
    mvBNode[nNumOfVert+iedge]=pBNode;
}
void CBoundaryTetra::deleteProgData()
{
    if(mnOrder==ElementOrder::First) {
        vector<CBoundaryNode*>().swap(mvEdgeBNode);
        vector<CBoundaryNode*>().swap(mvFaceBNode);
    } else {
        vector<CBoundaryNode*>().swap(mvFaceBNode);
    }
}
