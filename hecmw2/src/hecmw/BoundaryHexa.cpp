/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryHexa.cpp
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
#include <vector>
#include "BoundaryHexa.h"
#include "ElementProperty.h"
using namespace pmw;
CBoundaryHexa::CBoundaryHexa()
{
    uiint i;
    mvbMarkingEdge = new bool[NumberOfEdge::Hexa()];
    for(i=0; i < NumberOfEdge::Hexa(); i++) {
        mvbMarkingEdge[i]=false;
    };
    mvEdgeBNode.resize(NumberOfEdge::Hexa());
    mvEdgeNeibVol.resize(NumberOfEdge::Hexa());
    mvbMarkingFace = new bool[NumberOfFace::Hexa()];
    for(i=0; i < NumberOfFace::Hexa(); i++) {
        mvbMarkingFace[i]=false;
    };
    mvFaceBNode.resize(NumberOfFace::Hexa());
    mvFaceNeibVol.resize(NumberOfFace::Hexa());
    for(i=0; i < NumberOfFace::Hexa(); i++) {
        mvFaceBNode[i]= NULL;
    }
    mpVolBNode= NULL;
}
CBoundaryHexa::~CBoundaryHexa()
{
    ;
}
uiint CBoundaryHexa::getElemType()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    switch(mnOrder) {
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
    switch(mnOrder) {
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
    switch(mnOrder) {
    case(ElementOrder::First):
        mvBNode.resize(NumberOfNode::Hexa());
        break;
    case(ElementOrder::Second):
        mvBNode.resize(NumberOfNode::Hexa2());
        break;
    }
}
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
    return pairBNode;
}
uiint& CBoundaryHexa::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uiint ivert= getVertIndex(pairBNode.first);
    uiint jvert= getVertIndex(pairBNode.second);
    return pEdgeTree->getHexaEdgeIndex(ivert,jvert);
}
vector<CBoundaryNode*> CBoundaryHexa::getFaceCnvNodes(const uiint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    CBoundaryNode *pBNode;
    uiint  ivert, i, numOfVert;
    uiint *pvIndex;
    pvIndex= pFaceTree->getLocalNodeHexaFace(iface);
    numOfVert= pFaceTree->getHexaFaceNumOfVert(iface);
    vector<CBoundaryNode*> vFaceCnvNodes;
    vFaceCnvNodes.reserve(numOfVert);
    for(i=0; i < numOfVert; i++) {
        ivert= pvIndex[i];
        pBNode= mvBNode[ivert];
        vFaceCnvNodes.push_back(pBNode);
    }
    return vFaceCnvNodes;
}
uiint& CBoundaryHexa::getFaceID(vector<CBoundaryNode*>& vBNode)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    vuint vLocalVert;
    CBoundaryNode *pBNode;
    uiint ibnode, numOfBNode= vBNode.size();
    vLocalVert.reserve(4);
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= vBNode[ibnode];
        vLocalVert.push_back(getVertIndex(pBNode));
    };
    return pFaceTree->getHexaFaceIndex2(vLocalVert);
}
uiint* CBoundaryHexa::getLocalNode_Edge(const uiint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getHexaLocalNodeNum(iedge);
}
uiint* CBoundaryHexa::getLocalNode_Face(const uiint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    return pFaceTree->getLocalNodeHexaFace(iface);
}
void CBoundaryHexa::refine(uiint& countID, const vuint& vDOF)
{
    CBoundaryVolume *pProgVol;
    CElement   *pProgElem;
    CNode      *pNode;
    CBoundaryNode *pBNode;
    uiint iprog, numOfProg;
    uiint nProgPos;
    double progCubicVol;
    double coef;
    numOfProg= 8;
    mvProgVolume.reserve(numOfProg);
    for(iprog=0; iprog < numOfProg; iprog++) {
        pProgVol = new CBoundaryHexa;
        pProgVol->setID(countID);
        pProgVol->setOrder(mnOrder);
        countID++;
        mvProgVolume.push_back(pProgVol);
        nProgPos = dividHexa(iprog, pProgVol);
        pBNode= mvBNode[nProgPos];
        pNode = pBNode->getNode();
        pProgElem = mpElement->getProgElem_NodeID(pNode->getID());
        pProgVol->setElement(pProgElem);
        pProgVol->setElementID(pProgElem->getID());
        progCubicVol= pProgVol->calcVolume();
        coef = progCubicVol/mCubicVolume;
        distValue(pProgVol, coef, vDOF);
    };
}
double& CBoundaryHexa::calcVolume()
{
    CDiscreteVolume *pDiscrete= CDiscreteVolume::Instance();
    uiint* discreHexa;
    uiint  i,ii;
    CNode* vNode[4];
    mCubicVolume= 0.0;
    for(i=0; i< 6; i++) {
        discreHexa= pDiscrete->getHexaDiscrete(i);
        for(ii=0; ii< 4; ii++) {
            vNode[ii]= mvBNode[discreHexa[ii]]->getNode();
        };
        mCubicVolume += tetraVolume(vNode[0], vNode[1], vNode[2], vNode[3]);
    };
    return mCubicVolume;
}
void CBoundaryHexa::distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();
    double dAveVal, calcVal;
    uiint iedge, *pnEdgeVert;
    uiint nNumOfEdge=NumberOfEdge::Hexa();
    // 辺中心
    for(iedge=0; iedge < nNumOfEdge; iedge++) {
        pnEdgeVert = pEdgeTree->getHexaLocalNodeNum(iedge);
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
            mvEdgeBNode[iedge]->setValue(dof, mgLevel, calcVal);//-----------Level

            if(mgLevel!=nMaxMGLevel) {
                //mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);
                mvEdgeBNode[iedge]->setEntValue(dof, mgLevel+1, dAveVal);

                calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
                mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, calcVal);//--------Level+1
            }
        } else {
            if(mgLevel!=nMaxMGLevel) {
                //mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);
                mvEdgeBNode[iedge]->setEntValue(dof, mgLevel+1, dAveVal);

                calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
                mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, calcVal);//--------Level+1
            }
        }
    };

    CFaceTree *pFaceTree = CFaceTree::Instance();
    uiint iface, *pnFaceVert, ivert;
    uiint nNumOfFace = NumberOfFace::Hexa();
    // 面中心
    for(iface=0; iface < nNumOfFace; iface++) {
        pnFaceVert = pFaceTree->getLocalNodeHexaFace(iface);
        dAveVal=0.0;
        for(ivert=0; ivert < 4; ivert++) {
            //dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
            dAveVal += mvBNode[pnFaceVert[ivert]]->getEntValue(dof, mgLevel);
        };
        dAveVal *= 0.25;
        if(mgLevel!=nMaxMGLevel) {
            //mvFaceBNode[iface]->setValue(dof, mgLevel+1, dAveVal);
            mvFaceBNode[iface]->setEntValue(dof, mgLevel+1, dAveVal);

            calcVal= getCalcValue( mvFaceBNode[iface], dAveVal, dof, vPolandDOF, mPoland );
            mvFaceBNode[iface]->setValue(dof, mgLevel+1, calcVal);//-------Level+1
        }
    };
    // 頂点
    double dVertVal;
    dAveVal=0.0;
    uiint nNumOfVert=NumberOfVertex::Hexa();
    for(ivert=0; ivert < nNumOfVert; ivert++) {
        //dVertVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dVertVal = mvBNode[ivert]->getEntValue(dof, mgLevel);
        dAveVal += dVertVal;
        if(mgLevel!=nMaxMGLevel) {
            //mvBNode[ivert]->setValue(dof, mgLevel+1, dVertVal);
            mvBNode[ivert]->setEntValue(dof, mgLevel+1, dVertVal);

            calcVal= getCalcValue( mvBNode[ivert], dVertVal, dof, vPolandDOF, mPoland );
            mvBNode[ivert]->setValue(dof, mgLevel+1, calcVal);//-----Level+1
        }
    };
    // 要素中心
    dAveVal *= 0.125;
    if(mgLevel!=nMaxMGLevel) {
        //mpVolBNode->setValue(dof, mgLevel+1, dAveVal);
        mpVolBNode->setEntValue(dof, mgLevel+1, dAveVal);

        calcVal= getCalcValue( mpVolBNode, dAveVal, dof, vPolandDOF, mPoland );
        mpVolBNode->setValue(dof, mgLevel+1, calcVal);//----------Level+1
    }
}
void CBoundaryHexa::replaceEdgeBNode(const uiint& iedge)
{
    uiint nNumOfVert=NumberOfVertex::Hexa();
    CBoundaryNode *pBNode=mvEdgeBNode[iedge];
    mvBNode[nNumOfVert+iedge]=pBNode;
}
void CBoundaryHexa::deleteProgData()
{
    if(mnOrder==ElementOrder::First) {
        vector<CBoundaryNode*>().swap(mvEdgeBNode);
        vector<CBoundaryNode*>().swap(mvFaceBNode);
    } else {
        vector<CBoundaryNode*>().swap(mvFaceBNode);
    }
}
