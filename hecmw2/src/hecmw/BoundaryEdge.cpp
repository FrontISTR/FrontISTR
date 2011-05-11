
#include <vector>

#include "BoundaryParts.h"

//
//  BoundaryEdge.cpp
//
//
//
//              2010.04.13
//              k.Takeda
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

    // 1次-2次 Order
    switch(mnShapeType){
        case(ElementType::Beam):case(ElementType::Line):
            mnOrder = ElementOrder::First;
            break;
        case(ElementType::Beam2):case(ElementType::Line2):
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

// 辺Nodeを,辺BoundaryNodeにセット
// ----
// *mpElementは,CBounaryPartsが所有
// ----
void CBoundaryEdge::setupNode()
{
    CNode *pEdgeNode;

    pEdgeNode= mpElement->getEdgeInterNode(mnElemEdgeID);

    mpEdgeBNode->setNode(pEdgeNode);
}


// 辺の再分割
// ----
void CBoundaryEdge::refine(uiint& countID, const vuint& vDOF)
{
    mvProgBEdge.reserve(2);
    
    uiint ivert;
    for(ivert=0; ivert < 2; ivert++){

        CBoundaryEdge *pProgBEdge= new CBoundaryEdge;// <<<<<<<<<<<<<<<<< new
        
        mvProgBEdge.push_back(pProgBEdge);
        
        // 新BoundaryEdgeの属性
        pProgBEdge->resizeBNode(2);
        pProgBEdge->setID(countID);
        countID++;

        // タイプをセット
        if(mnOrder==ElementOrder::First) pProgBEdge->setBEdgeShape(ElementType::Beam);
        if(mnOrder==ElementOrder::Second)pProgBEdge->setBEdgeShape(ElementType::Beam2);

        // 要素のNodeから要素頂点番号を取得 -> 子要素を取得
        CNode *pNode= mvBNode[ivert]->getNode();
        uiint elemVert= mpElement->getLocalVertNum(pNode->getID());
        CElement *pProgElem= mpElement->getProgElem(elemVert);

        pProgBEdge->setElement(pProgElem);
        pProgBEdge->setElementID(pProgElem->getID());

        // 2点のNodeから要素の辺番号を取得
        CNode *pEdgeNode= mpEdgeBNode->getNode();
        uiint elemEdge= pProgElem->getEdgeIndex(pNode, pEdgeNode);//2010.07.22
        pProgBEdge->setElementEdgeID(elemEdge);
        
        
        // BNodeをセット
        pProgBEdge->setBNode(0, mvBNode[ivert]);
        pProgBEdge->setBNode(1, mpEdgeBNode);
        
        // 線分比による境界値の分配
        double coef= pProgBEdge->calcLength() / mLength;
        
        ////debug
        //cout << "BoundaryEdge::refine, mLength= " << mLength << ", coef= " << coef << endl;

        double progValue;
        uiint idof, dof;
        for(idof=0; idof < vDOF.size(); idof++){
            dof = vDOF[idof];
            progValue = mmValue[dof]*coef;

            pProgBEdge->setBndValue(dof,progValue);
        };
    };
}


// 自身の線分の距離
// ----
double& CBoundaryEdge::calcLength()
{
    CNode *pNode0, *pNode1;
    pNode0= mvBNode[0]->getNode();
    pNode1= mvBNode[1]->getNode();
    double x, y, z;
    
    x = pNode0->getX() - pNode1->getX();
    y = pNode0->getY() - pNode1->getY();
    z = pNode0->getZ() - pNode1->getZ();

    ////debug
    ////cout << "BoundaryEdge::calcLength, x= " << x << ", y= " << y << ", z= " << z << endl;
    //cout << "BoundaryEdge::calcLength, Node0->ID= " << pNode0->getID() << ", Node1->ID= " << pNode1->getID() << endl;
    //cout << "BoundaryEdge::calcLength, Node0->X= " << pNode0->getX() << ", Node1->X= " << pNode1->getX() << endl;
    //cout << "BoundaryEdge::calcLength, Node0->Y= " << pNode0->getY() << ", Node1->Y= " << pNode1->getY() << endl;
    //cout << "BoundaryEdge::calcLength, Node0->Z= " << pNode0->getZ() << ", Node1->Z= " << pNode1->getZ() << endl;

    mLength = sqrt(x*x + y*y + z*z);
    
    return mLength;
}


// 2次要素の場合、辺BNodeをmvBNodeに移設.
//
void CBoundaryEdge::replaceEdgeBNode()
{
    if(mnOrder==ElementOrder::Second){
        mvBNode.resize(3);
        mvBNode[2] = mpEdgeBNode;
    }

    //cout << "BoundaryEdge::replaceEdgeBNode, mvBNode.size " << mvBNode.size() << endl;
}


// 上位GridのBNodeへディレクレ値をセット
//
void CBoundaryEdge::distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel)
{
    double dAveVal(0.0), dVal;

    dVal = mvBNode[0]->getValue(dof, mgLevel);
    if(mgLevel!=nMaxMGLevel){
      mvBNode[0]->setValue(dof, mgLevel+1, dVal);//上位グリッドへ下位の値をセット:頂点
    }
    dAveVal += dVal;

    dVal = mvBNode[1]->getValue(dof, mgLevel);
    if(mgLevel!=nMaxMGLevel){
      mvBNode[1]->setValue(dof, mgLevel+1, dVal);//上位グリッドへ下位の値をセット:頂点
    }
    dAveVal += dVal;

    // 辺中心値
    dAveVal *= 0.5;
    if(mnOrder==ElementOrder::Second){
        //2次要素:カレント・グリッドと上位グリッド
        mpEdgeBNode->setValue(dof, mgLevel,   dAveVal);
        if(mgLevel!=nMaxMGLevel)
          mpEdgeBNode->setValue(dof, mgLevel+1, dAveVal);//辺中央のBNodeは、上位Gridへディレクレ値をセット
    }else{
        //1次要素:上位グリッド
        if(mgLevel!=nMaxMGLevel)
          mpEdgeBNode->setValue(dof, mgLevel+1, dAveVal);//辺中央のBNodeは、上位Gridへディレクレ値をセット
    }

////    if(mgLevel > 0){
////        dAveVal *= 0.5;
////
////        //辺BNodeへのディレクレ値
////        if(mnOrder==ElementOrder::Second){
////            //2次要素:カレント・グリッドと上位グリッド
////            mpEdgeBNode->setValue(dof, mgLevel,   dAveVal);
////            if(mgLevel!=nMaxMGLevel)
////              mpEdgeBNode->setValue(dof, mgLevel+1, dAveVal);//辺中央のBNodeは、上位Gridへディレクレ値をセット
////        }else{
////            //1次要素:上位グリッド
////            if(mgLevel!=nMaxMGLevel)
////              mpEdgeBNode->setValue(dof, mgLevel+1, dAveVal);//辺中央のBNodeは、上位Gridへディレクレ値をセット
////        }
////    }
////    if(mgLevel==0){
////        if(mnOrder==ElementOrder::Second){
////            //2次要素:カレント・グリッドと上位グリッド
////            mpEdgeBNode->setValue(dof, mgLevel,   mmValue[dof]);
////            if(mgLevel!=nMaxMGLevel)
////              mpEdgeBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
////        }else{
////            //1次要素:上位グリッド
////            if(mgLevel!=nMaxMGLevel)
////              mpEdgeBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
////        }
////    }
}




