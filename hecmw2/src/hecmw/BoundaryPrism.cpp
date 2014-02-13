/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryPrism.cpp
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
#include "BoundaryPrism.h"
#include "ElementProperty.h"
using namespace pmw;
CBoundaryPrism::CBoundaryPrism()
{
    uiint i;
    mvbMarkingEdge = new bool[NumberOfEdge::Prism()];
    for(i=0; i < NumberOfEdge::Prism(); i++) {
        mvbMarkingEdge[i]=false;
    };
    mvEdgeBNode.resize(NumberOfEdge::Prism());
    mvEdgeNeibVol.resize(NumberOfEdge::Prism());
    mvbMarkingFace = new bool[NumberOfFace::Prism()];
    for(i=0; i < NumberOfFace::Prism(); i++) {
        mvbMarkingFace[i]=false;
    };
    mvFaceBNode.resize(NumberOfFace::Prism());
    mvFaceNeibVol.resize(NumberOfFace::Prism());
    for(i=0; i < NumberOfFace::Prism(); i++) {
        mvFaceBNode[i]= NULL;
    }
    mpVolBNode= NULL;
}
CBoundaryPrism::~CBoundaryPrism()
{
    ;
}
uiint CBoundaryPrism::getElemType()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    switch(mnOrder) {
    case(ElementOrder::First):
        return ElementType::Prism;
    case(ElementOrder::Second):
        return ElementType::Prism2;
    default:
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryPrism::getElemType, invalid mnOrder");
        return ElementType::Limit;
    }
}
uiint CBoundaryPrism::getNumOfEdge()
{
    return NumberOfEdge::Prism();
}
uiint CBoundaryPrism::getNumOfFace()
{
    return NumberOfFace::Prism();
}
uiint CBoundaryPrism::getNumOfNode()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    switch(mnOrder) {
    case(ElementOrder::First):
        return NumberOfNode::Prism();
    case(ElementOrder::Second):
        return NumberOfNode::Prism2();
    default:
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryPrism::getNumOfNode, invalid mnOrder");
        return NumberOfNode::Default();
    }
}
uiint CBoundaryPrism::getNumOfVert()
{
    return NumberOfVertex::Prism();;
}
void CBoundaryPrism::setOrder(const uiint& order)
{
    mnOrder = order;
    switch(mnOrder) {
    case(ElementOrder::First):
        mvBNode.resize(NumberOfNode::Prism());
        break;
    case(ElementOrder::Second):
        mvBNode.resize(NumberOfNode::Prism2());
        break;
    }
}
PairBNode CBoundaryPrism::getPairBNode(const uiint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uiint *pLocalNum;
    pLocalNum= pEdgeTree->getPrismLocalNodeNum(iedge);
    uiint index1st = pLocalNum[0];
    uiint index2nd = pLocalNum[1];
    PairBNode pairBNode;
    pairBNode.first = mvBNode[index1st];
    pairBNode.second= mvBNode[index2nd];
    return pairBNode;
}
uiint& CBoundaryPrism::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uiint ivert= getVertIndex(pairBNode.first);
    uiint jvert= getVertIndex(pairBNode.second);
    return pEdgeTree->getPrismEdgeIndex(ivert, jvert);
}
vector<CBoundaryNode*> CBoundaryPrism::getFaceCnvNodes(const uiint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    CBoundaryNode *pBNode;
    uiint  ivert, i, numOfVert;
    uiint *pvIndex;
    pvIndex= pFaceTree->getLocalNodePrismFace(iface);
    numOfVert= pFaceTree->getPrismFaceNumOfVert(iface);
    vector<CBoundaryNode*> vFaceCnvNodes;
    vFaceCnvNodes.reserve(numOfVert);
    for(i=0; i < numOfVert; i++) {
        ivert= pvIndex[i];
        pBNode= mvBNode[ivert];
        vFaceCnvNodes.push_back(pBNode);
    }
    return vFaceCnvNodes;
}
uiint& CBoundaryPrism::getFaceID(vector<CBoundaryNode*>& vBNode)
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
    return pFaceTree->getPrismFaceIndex2(vLocalVert);
}
uiint* CBoundaryPrism::getLocalNode_Edge(const uiint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getPrismLocalNodeNum(iedge);
}
uiint* CBoundaryPrism::getLocalNode_Face(const uiint& iface)
{
    CFaceTree *pFaceTree= CFaceTree::Instance();
    return pFaceTree->getLocalNodePrismFace(iface);
}
void CBoundaryPrism::refine(uiint& countID, const vuint& vDOF)
{
    CBoundaryVolume *pProgVol;
    CElement   *pProgElem;
    CNode      *pNode;
    CBoundaryNode *pBNode;
    uiint iprog, numOfProg;
    uiint nProgPos;
    double progCubicVol;
    double coef;
    numOfProg= 6;
    mvProgVolume.reserve(numOfProg);
    for(iprog=0; iprog < numOfProg; iprog++) {
        pProgVol = new CBoundaryHexa;
        pProgVol->setID(countID);
        pProgVol->setOrder(mnOrder);
        countID++;
        mvProgVolume.push_back(pProgVol);
        nProgPos = dividPrism(iprog, pProgVol);
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
double& CBoundaryPrism::calcVolume()
{
    CDiscreteVolume *pDiscrete= CDiscreteVolume::Instance();
    uiint* discrePrism;
    uiint  i,ii;
    CNode* vNode[4];
    mCubicVolume= 0.0;
    for(i=0; i< 3; i++) {
        discrePrism= pDiscrete->getPrismDiscrete(i);
        for(ii=0; ii< 4; ii++) {
            vNode[ii]= mvBNode[discrePrism[ii]]->getNode();
        };
        mCubicVolume += tetraVolume(vNode[0], vNode[1], vNode[2], vNode[3]);
    };
    return mCubicVolume;
}
void CBoundaryPrism::distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();
    double dAveVal, calcVal;
    uiint iedge, *pnEdgeVert;
    uiint nNumOfEdge=NumberOfEdge::Prism();
    // 辺中心
    for(iedge=0; iedge < nNumOfEdge; iedge++) {
        pnEdgeVert = pEdgeTree->getPrismLocalNodeNum(iedge);
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
            mvEdgeBNode[iedge]->setValue(dof, mgLevel, calcVal);//----------Level
        }

        if(mgLevel!=nMaxMGLevel) {
            //mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);
            mvEdgeBNode[iedge]->setEntValue(dof, mgLevel+1, dAveVal);

            calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
            mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, calcVal);//--------Level+1
        }

    };
    // 面中心
    CFaceTree *pFaceTree = CFaceTree::Instance();
    uiint iface, *pnFaceVert, ivert;
    uiint nNumOfFace=NumberOfFace::Prism();
    for(iface=0; iface < nNumOfFace; iface++) {
        pnFaceVert = pFaceTree->getLocalNodePrismFace(iface);
        dAveVal=0.0;
        if(iface > 1) {
            for(ivert=0; ivert < 4; ivert++) {
                //dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
                dAveVal += mvBNode[pnFaceVert[ivert]]->getEntValue(dof, mgLevel);
            };
            dAveVal *= 0.25;
        } else {
            for(ivert=0; ivert < 3; ivert++) {
                //dAveVal += mvBNode[pnFaceVert[ivert]]->getValue(dof, mgLevel);
                dAveVal += mvBNode[pnFaceVert[ivert]]->getEntValue(dof, mgLevel);
            };
            dAveVal /= 3.0;
        }

        if(mgLevel!=nMaxMGLevel) {
            //mvFaceBNode[iface]->setValue(dof, mgLevel+1, dAveVal);
            mvFaceBNode[iface]->setEntValue(dof, mgLevel+1, dAveVal);

            calcVal= getCalcValue( mvFaceBNode[iface], dAveVal, dof, vPolandDOF, mPoland );
            mvFaceBNode[iface]->setValue(dof, mgLevel+1, calcVal);//--------Level+1
        }
    };
    // 頂点
    double dVertVal;
    dAveVal=0.0;
    uiint nNumOfVert=NumberOfVertex::Prism();
    for(ivert=0; ivert < nNumOfVert; ivert++) {
        //dVertVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dVertVal = mvBNode[ivert]->getEntValue(dof, mgLevel);
        dAveVal += dVertVal;
        if(mgLevel!=nMaxMGLevel) {
            //mvBNode[ivert]->setValue(dof, mgLevel+1, dVertVal);
            mvBNode[ivert]->setEntValue(dof, mgLevel+1, dVertVal);

            calcVal= getCalcValue( mvBNode[ivert], dVertVal, dof, vPolandDOF, mPoland );
            mvBNode[ivert]->setValue(dof, mgLevel+1, calcVal);//-----------Level+1
        }
    };
    // 要素中心
    dAveVal /= 6.0;
    if(mgLevel!=nMaxMGLevel) {
        //mpVolBNode->setValue(dof, mgLevel+1, dAveVal);
        mpVolBNode->setEntValue(dof, mgLevel+1, dAveVal);

        calcVal= getCalcValue( mpVolBNode, dAveVal, dof, vPolandDOF, mPoland );
        mpVolBNode->setValue(dof, mgLevel+1, calcVal);//---------------Level+1
    }
}
void CBoundaryPrism::replaceEdgeBNode(const uiint& iedge)
{
    uiint nNumOfVert=NumberOfVertex::Prism();
    CBoundaryNode *pBNode=mvEdgeBNode[iedge];
    mvBNode[nNumOfVert+iedge]=pBNode;
}
void CBoundaryPrism::deleteProgData()
{
    if(mnOrder==ElementOrder::First) {
        vector<CBoundaryNode*>().swap(mvEdgeBNode);
        vector<CBoundaryNode*>().swap(mvFaceBNode);
    } else {
        vector<CBoundaryNode*>().swap(mvFaceBNode);
    }
}
