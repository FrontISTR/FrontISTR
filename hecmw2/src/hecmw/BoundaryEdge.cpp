/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryEdge.cpp
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
#include "BoundaryParts.h"
#include "BoundaryEdge.h"
#include "ElementProperty.h"
using namespace pmw;
CBoundaryEdge::CBoundaryEdge()
{
    mpEdgeBNode = NULL;
}
CBoundaryEdge::~CBoundaryEdge()
{
    ;
}
void CBoundaryEdge::setBEdgeShape(const uiint& elemType)
{
    mnShapeType = elemType;
    switch(mnShapeType) {
    case(ElementType::Beam):
    case(ElementType::Line):
        mnOrder = ElementOrder::First;
        break;
    case(ElementType::Beam2):
    case(ElementType::Line2):
        mnOrder = ElementOrder::Second;
        break;
    }
}
uiint CBoundaryEdge::getNumOfVert()
{
    return 2;
}
void CBoundaryEdge::setEdgeBNode(CBoundaryNode* pBNode)
{
    mpEdgeBNode= pBNode;
}
PairBNode CBoundaryEdge::getPairBNode()
{
    PairBNode pairBNode;
    pairBNode.first = mvBNode[0];
    pairBNode.second= mvBNode[1];
    return pairBNode;
}
void CBoundaryEdge::setupNode()
{
    CNode *pEdgeNode;
    pEdgeNode= mpElement->getEdgeInterNode(mnElemEdgeID);
    mpEdgeBNode->setNode(pEdgeNode);
}
void CBoundaryEdge::refine(uiint& countID, const vuint& vDOF)
{
    mvProgBEdge.reserve(2);
    uiint ivert;
    for(ivert=0; ivert < 2; ivert++) {
        CBoundaryEdge *pProgBEdge= new CBoundaryEdge;
        mvProgBEdge.push_back(pProgBEdge);
        pProgBEdge->resizeBNode(2);
        pProgBEdge->setID(countID);
        countID++;
        if(mnOrder==ElementOrder::First) pProgBEdge->setBEdgeShape(ElementType::Beam);
        if(mnOrder==ElementOrder::Second)pProgBEdge->setBEdgeShape(ElementType::Beam2);
        CNode *pNode= mvBNode[ivert]->getNode();
        uiint elemVert= mpElement->getLocalVertNum(pNode->getID());
        CElement *pProgElem= mpElement->getProgElem(elemVert);
        pProgBEdge->setElement(pProgElem);
        pProgBEdge->setElementID(pProgElem->getID());
        CNode *pEdgeNode= mpEdgeBNode->getNode();
        uiint elemEdge= pProgElem->getEdgeIndex(pNode, pEdgeNode);
        pProgBEdge->setElementEdgeID(elemEdge);
        pProgBEdge->setBNode(0, mvBNode[ivert]);
        pProgBEdge->setBNode(1, mpEdgeBNode);
        double coef= pProgBEdge->calcLength() / mLength;
        double progValue;
        uiint idof, dof;
        for(idof=0; idof < vDOF.size(); idof++) {
            dof = vDOF[idof];
            progValue = mmValue[dof]*coef;
            pProgBEdge->setBndValue(dof,progValue);
        };
    };
}
double& CBoundaryEdge::calcLength()
{
    CNode *pNode0, *pNode1;
    pNode0= mvBNode[0]->getNode();
    pNode1= mvBNode[1]->getNode();
    double x, y, z;
    x = pNode0->getX() - pNode1->getX();
    y = pNode0->getY() - pNode1->getY();
    z = pNode0->getZ() - pNode1->getZ();
    mLength = sqrt(x*x + y*y + z*z);
    return mLength;
}
void CBoundaryEdge::replaceEdgeBNode()
{
    if(mnOrder==ElementOrder::Second) {
        mvBNode.resize(3);
        mvBNode[2] = mpEdgeBNode;
    }
}
void CBoundaryEdge::distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland)
{
    double dAveVal(0.0), dVal, calcVal;
    //--
    // 端点0
    //--
    //dVal = mvBNode[0]->getValue(dof, mgLevel);
    dVal = mvBNode[0]->getEntValue(dof, mgLevel);//カレントレベル

    if(mgLevel!=nMaxMGLevel) {
        //mvBNode[0]->setValue(dof, mgLevel+1, dVal);
        mvBNode[0]->setEntValue(dof, mgLevel+1, dVal);//上位レベル

        calcVal= getCalcValue( mvBNode[0], dVal, dof, vPolandDOF, mPoland );
        mvBNode[0]->setValue(dof, mgLevel+1, calcVal);//-----Level+1: 数式処理あり=>数式計算値，数式なし=>基礎データ
    }
    //--
    // 端点1
    //--
    dAveVal += dVal;
    //dVal = mvBNode[1]->getValue(dof, mgLevel);
    dVal = mvBNode[1]->getEntValue(dof, mgLevel);//カレントレベル

    if(mgLevel!=nMaxMGLevel) {
        //mvBNode[1]->setValue(dof, mgLevel+1, dVal);
        mvBNode[1]->setEntValue(dof, mgLevel+1, dVal);//上位レベル

        calcVal= getCalcValue( mvBNode[1], dVal, dof, vPolandDOF, mPoland );
        mvBNode[1]->setValue(dof, mgLevel+1, calcVal);//-----Level+1: 数式処理あり=>数式計算値，数式なし=>基礎データ
    }
    //--
    // 辺中心
    //--
    dAveVal += dVal;
    dAveVal *= 0.5;
    if(mnOrder==ElementOrder::Second) {
        //mpEdgeBNode->setValue(dof, mgLevel,   dAveVal);
        mpEdgeBNode->setEntValue(dof, mgLevel,   dAveVal);

        calcVal= getCalcValue( mpEdgeBNode, dAveVal, dof, vPolandDOF, mPoland );
        mpEdgeBNode->setValue(dof, mgLevel, calcVal);//-----Level: 数式処理あり=>数式計算値，数式なし=>基礎データ

        if(mgLevel!=nMaxMGLevel) {
            //mpEdgeBNode->setValue(dof, mgLevel+1, dAveVal);
            mpEdgeBNode->setEntValue(dof, mgLevel+1, dAveVal);

            calcVal= getCalcValue( mpEdgeBNode, dAveVal, dof, vPolandDOF, mPoland );
            mpEdgeBNode->setValue(dof, mgLevel+1, calcVal);//-----Level+1: 数式処理あり=>数式計算値，数式なし=>基礎データ
        }
    } else {
        if(mgLevel!=nMaxMGLevel) {
            //mpEdgeBNode->setValue(dof, mgLevel+1, dAveVal);
            mpEdgeBNode->setEntValue(dof, mgLevel+1, dAveVal);

            calcVal= getCalcValue( mpEdgeBNode, dAveVal, dof, vPolandDOF, mPoland );
            mpEdgeBNode->setValue(dof, mgLevel+1, calcVal);//-----Level+1: 数式処理あり=>数式計算値，数式なし=>基礎データ
        }
    }
}
