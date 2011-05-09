
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
using namespace pmw;

CBoundaryEdge::CBoundaryEdge()
{
    ;
}

CBoundaryEdge::~CBoundaryEdge()
{
    ;
}

void CBoundaryEdge::setEdgeBNode(CBoundaryNode* pBNode)
{
    mpEdgeBNode= pBNode;
}

PairBNode& CBoundaryEdge::getPairBNode()
{
    mPairBNode.first = mvBNode[0];
    mPairBNode.second= mvBNode[1];

    return mPairBNode;
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
void CBoundaryEdge::refine(uint& countID, const vuint& vDOF)
{
    mvProgBEdge.reserve(2);
    
    uint ivert;
    for(ivert=0; ivert < 2; ivert++){

        CBoundaryEdge *pProgBEdge= new CBoundaryEdge;// <<<<<<<<<<<<<<<<< new
        
        mvProgBEdge.push_back(pProgBEdge);
        
        // 新BoundaryEdgeの属性
        pProgBEdge->resizeBNode(2);
        pProgBEdge->setID(countID);
        countID++;
        pProgBEdge->setBEdgeShape(ElementType::Beam);// タイプをセット

        // 要素のNodeから要素頂点番号を取得 -> 子要素を取得
        CNode *pNode= mvBNode[ivert]->getNode();
        uint elemVert= mpElement->getLocalVertNum(pNode->getID());
        CElement *pProgElem= mpElement->getProgElem(elemVert);

        pProgBEdge->setElement(pProgElem);
        pProgBEdge->setElementID(pProgElem->getID());

        // 2点のNodeから要素の辺番号を取得
        CNode *pEdgeNode= mpEdgeBNode->getNode();
        uint elemEdge= pProgElem->getEdgeIndex(pNode, pEdgeNode);//2010.07.22 
        pProgBEdge->setElementEdgeID(elemEdge);
        
        
        // BNodeをセット
        pProgBEdge->setBNode(0, mvBNode[ivert]);
        pProgBEdge->setBNode(1, mpEdgeBNode);
        
        // 線分比による境界値の分配
        double coef= pProgBEdge->calcLength() / mLength;
        
        ////debug
        //cout << "BoundaryEdge::refine, mLength= " << mLength << ", coef= " << coef << endl;

        double progValue;
        uint idof, dof;
        for(idof=0; idof < vDOF.size(); idof++){
            dof = vDOF[idof];
            progValue = mmValue[dof]*coef;

            //pProgBEdge->addDOF(dof);
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


// 上位GridのBNodeへディレクレ値をセット
//
void CBoundaryEdge::distDirichletVal(const uint& dof, const uint& mgLevel)
{
    double dAveVal(0.0), dVal;
    
    dVal = mvBNode[0]->getValue(dof, mgLevel);
    mvBNode[0]->setValue(dof, mgLevel+1, dVal);//上位グリッドへ下位の値をセット:頂点
    dAveVal += dVal;

    dVal = mvBNode[1]->getValue(dof, mgLevel);
    mvBNode[1]->setValue(dof, mgLevel+1, dVal);//上位グリッドへ下位の値をセット:頂点
    dAveVal += dVal;

    if(mgLevel > 0){
        dAveVal *= 0.5;
        mpEdgeBNode->setValue(dof, mgLevel+1, dAveVal);//辺中央のBNodeは、上位Gridへディレクレ値をセット
    }
    if(mgLevel==0){
        mpEdgeBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
    }
}




