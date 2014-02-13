/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryFace.cpp
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
#include "EdgeTree.h"
#include "BoundaryNode.h"
#include "BoundaryParts.h"
#include "BoundaryFace.h"
#include "ElementType.h"
#include "HEC_MPI.h"
using namespace pmw;
CBoundaryFace::CBoundaryFace()
{
    mArea = 0.0;
    mpFaceBNode= NULL;
}
CBoundaryFace::~CBoundaryFace()
{
    delete []mvbMarkingEdge;
}
uiint CBoundaryFace::getNumOfEdge()
{
    switch(mnShapeType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        return NumberOfEdge::Quad();
        break;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        return NumberOfEdge::Triangle();
        break;
    default:
        return NumberOfEdge::Default();
        break;
    }
}
void CBoundaryFace::setBFaceShape(const uiint& elemType)
{
    Utility::CLogger* pLogger = Utility::CLogger::Instance();
    mnShapeType= elemType;
    uiint ivert;
    switch(mnShapeType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        mvEdgeNeibFace.resize(4);
        mvbMarkingEdge = new bool[4];
        for(ivert=0; ivert < 4; ivert++)
            mvbMarkingEdge[ivert] = false;
        mvEdgeBNode.resize(4);
        break;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        mvEdgeNeibFace.resize(3);
        mvbMarkingEdge = new bool[3];
        for(ivert=0; ivert < 3; ivert++)
            mvbMarkingEdge[ivert] = false;
        mvEdgeBNode.resize(3);
        break;
    default:
        pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::setBFaceShape");
        break;
    }
    switch(mnShapeType) {
    case(ElementType::Quad):
    case(ElementType::Triangle):
        mnOrder = ElementOrder::First;
        break;
    case(ElementType::Quad2):
    case(ElementType::Triangle2):
        mnOrder = ElementOrder::Second;
        break;
    }
}
uiint CBoundaryFace::getNumOfVert()
{
    switch(mnShapeType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        return 4;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        return 3;
    default:
        return 0;
    }
}
void CBoundaryFace::markingEdge(const uiint& iedge)
{
    mvbMarkingEdge[iedge]= true;
}
void CBoundaryFace::setEdgeNeibFace(const uiint& iedge, const uiint& neibFaceID)
{
    mvEdgeNeibFace[iedge]= neibFaceID;
}
void CBoundaryFace::setEdgeBNode(const uiint& iedge, CBoundaryNode* pEdgeBNode)
{
    mvEdgeBNode[iedge]= pEdgeBNode;
}
void CBoundaryFace::setFaceBNode(CBoundaryNode* pFaceBNode)
{
    mpFaceBNode= pFaceBNode;
}
PairBNode CBoundaryFace::getPairBNode(const uiint& iedge)
{
    Utility::CLogger* pLogger = Utility::CLogger::Instance();
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uiint *pLocalNum;
    switch(mnShapeType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        pLocalNum= pEdgeTree->getQuadLocalNodeNum(iedge);
        break;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        pLocalNum= pEdgeTree->getTriangleLocalNodeNum(iedge);
        break;
    default:
        pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::getPairBNode");
        break;
    }
    uiint index1st = pLocalNum[0];
    uiint index2nd = pLocalNum[1];
    PairBNode pairBNode;
    pairBNode.first = mvBNode[index1st];
    pairBNode.second= mvBNode[index2nd];
    return pairBNode;
}
uiint& CBoundaryFace::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uiint ivert= getVertIndex(pairBNode.first);
    uiint jvert= getVertIndex(pairBNode.second);
    switch(mnShapeType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        return pEdgeTree->getQuadEdgeIndex(ivert,jvert);
        break;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        return pEdgeTree->getTriangleEdgeIndex(ivert, jvert);
        break;
    default:
        return pEdgeTree->getDisagTypeEdgeIndex(ivert, jvert);
        break;
    }
}
void CBoundaryFace::setupNode_Edge()
{
    uiint nNumOfEdge = getNumOfEdge();
    CNode *pEdgeNode, *pNode0, *pNode1;
    uiint elemEdgeID;
    CBoundaryNode *pEdgeBNode, *pBNode0, *pBNode1;
    uiint iedge;
    for(iedge=0; iedge < nNumOfEdge; iedge++) {
        pEdgeBNode= mvEdgeBNode[iedge];
        if(iedge != nNumOfEdge-1) {
            pBNode0= mvBNode[iedge];
            pBNode1= mvBNode[iedge+1];
        } else {
            pBNode0= mvBNode[iedge];
            pBNode1= mvBNode[0];
        }
        pNode0= pBNode0->getNode();
        pNode1= pBNode1->getNode();
        elemEdgeID= mpElement->getEdgeIndex(pNode0, pNode1);
        pEdgeNode = mpElement->getEdgeInterNode(elemEdgeID);
        pEdgeBNode->setNode(pEdgeNode);
    };
}
void CBoundaryFace::setupNode_Face()
{
    CNode *pFaceNode;
    pFaceNode= mpElement->getFaceNode(mnElemFaceID);
    mpFaceBNode->setNode(pFaceNode);
}

void CBoundaryFace::refine(uiint& countID, const vuint& vDOF)
{
    Utility::CLogger* pLogger = Utility::CLogger::Instance();
    uiint numOfProg;
    uiint iprog;
    CBoundaryFace *pProgBFace;
    double progArea;
    double progValue;
    double coef;
    uiint idof, dof;
    CElement  *pProgElem;
    uiint     nProgEntity;
    uiint        ipos;
    CNode      *pNode;
    vector<CNode*>  vNode;
    CBoundaryNode *pBNode;
    uiint          ibnode;

    switch(mnShapeType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        numOfProg= 4;
        mvProgBFace.reserve(numOfProg);
        for(iprog=0; iprog < numOfProg; iprog++) {
            pProgBFace = new CBoundaryFace;//--------------------------------生成
            pProgBFace->setID(countID);
            countID++;
            if(mnOrder==ElementOrder::First ) pProgBFace->setBFaceShape(ElementType::Quad);
            if(mnOrder==ElementOrder::Second) pProgBFace->setBFaceShape(ElementType::Quad2);
            pProgBFace->resizeBNode(4);
            switch(iprog) {
            case(0):
                pProgBFace->setBNode(0, mvBNode[0]);
                pProgBFace->setBNode(1, mvEdgeBNode[0]);
                pProgBFace->setBNode(2, mpFaceBNode);
                pProgBFace->setBNode(3, mvEdgeBNode[3]);
                progArea= pProgBFace->calcArea();
                ;
                break;
            case(1):
                pProgBFace->setBNode(0, mvEdgeBNode[0]);
                pProgBFace->setBNode(1, mvBNode[1]);
                pProgBFace->setBNode(2, mvEdgeBNode[1]);
                pProgBFace->setBNode(3, mpFaceBNode);
                progArea= pProgBFace->calcArea();
                break;
            case(2):
                pProgBFace->setBNode(0, mvEdgeBNode[1]);
                pProgBFace->setBNode(1, mvBNode[2]);
                pProgBFace->setBNode(2, mvEdgeBNode[2]);
                pProgBFace->setBNode(3, mpFaceBNode);
                progArea= pProgBFace->calcArea();
                break;
            case(3):
                pProgBFace->setBNode(0, mvEdgeBNode[2]);
                pProgBFace->setBNode(1, mvBNode[3]);
                pProgBFace->setBNode(2, mvEdgeBNode[3]);
                pProgBFace->setBNode(3, mpFaceBNode);
                progArea= pProgBFace->calcArea();
                break;
            }
            mvProgBFace.push_back(pProgBFace);
            coef = progArea / mArea;
            for(idof=0; idof < vDOF.size(); idof++) {
                dof = vDOF[idof];
                progValue = mmValue[dof]*coef;//----------------基礎境界値計算(面が持っている値):Neumannの元
                pProgBFace->setBndValue(dof,progValue);
            };
            pNode = mvBNode[iprog]->getNode();
            ipos  = mpElement->getLocalVertNum(pNode->getID());
            pProgElem = mpElement->getProgElem(ipos);
            pProgBFace->setElement(pProgElem);
            pProgBFace->setElementID(pProgElem->getID());
            vNode.clear();
            for(ibnode=0; ibnode < 3; ibnode++) {
                pBNode = pProgBFace->getBNode(ibnode);
                pNode = pBNode->getNode();
                vNode.push_back(pNode);
            };
            nProgEntity = pProgElem->getFaceIndex(vNode[0], vNode[1], vNode[2]);
            pProgBFace->setElementFaceID(nProgEntity);
        };
        break;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        numOfProg= 3;
        mvProgBFace.reserve(numOfProg);
        for(iprog=0; iprog < numOfProg; iprog++) {
            pProgBFace = new CBoundaryFace;//-----------------------------------------------生成
            pProgBFace->setID(countID);
            countID++;
            if(mnOrder==ElementOrder::First ) pProgBFace->setBFaceShape(ElementType::Quad);
            if(mnOrder==ElementOrder::Second) pProgBFace->setBFaceShape(ElementType::Quad2);
            pProgBFace->resizeBNode(4);
            switch(iprog) {
            case(0):
                pProgBFace->setBNode(0, mvBNode[0]);
                pProgBFace->setBNode(1, mvEdgeBNode[0]);
                pProgBFace->setBNode(2, mpFaceBNode);
                pProgBFace->setBNode(3, mvEdgeBNode[2]);
                progArea= pProgBFace->calcArea();
                break;
            case(1):
                pProgBFace->setBNode(0, mvEdgeBNode[0]);
                pProgBFace->setBNode(1, mvBNode[1]);
                pProgBFace->setBNode(2, mvEdgeBNode[1]);
                pProgBFace->setBNode(3, mpFaceBNode);
                progArea= pProgBFace->calcArea();
                break;
            case(2):
                pProgBFace->setBNode(0, mvEdgeBNode[1]);
                pProgBFace->setBNode(1, mvBNode[2]);
                pProgBFace->setBNode(2, mvEdgeBNode[2]);
                pProgBFace->setBNode(3, mpFaceBNode);
                progArea= pProgBFace->calcArea();
                break;
            }
            mvProgBFace.push_back(pProgBFace);
            coef = progArea / mArea;
            for(idof=0; idof < vDOF.size(); idof++) {
                dof = vDOF[idof];
                progValue = mmValue[dof]*coef;//-------------------基礎境界値計算(面が持っている値):Neumannの元
                pProgBFace->setBndValue(dof,progValue);
            };
            pNode = mvBNode[iprog]->getNode();
            ipos  = mpElement->getLocalVertNum(pNode->getID());
            pProgElem = mpElement->getProgElem(ipos);
            pProgBFace->setElement(pProgElem);
            pProgBFace->setElementID(pProgElem->getID());
            vNode.clear();
            for(ibnode=0; ibnode < 3; ibnode++) {
                pBNode = pProgBFace->getBNode(ibnode);
                pNode = pBNode->getNode();
                vNode.push_back(pNode);
            };
            nProgEntity = pProgElem->getFaceIndex(vNode[0], vNode[1], vNode[2]);
            pProgBFace->setElementFaceID(nProgEntity);
        };
        break;
    default:
        pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::refine");
        break;
    }
}
double CBoundaryFace::triArea(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    double x, y, z, x1, y1, z1, x2, y2, z2;
    x1 = pNode1->getX() - pNode0->getX();
    y1 = pNode1->getY() - pNode0->getY();
    z1 = pNode1->getZ() - pNode0->getZ();
    x2 = pNode2->getX() - pNode0->getX();
    y2 = pNode2->getY() - pNode0->getY();
    z2 = pNode2->getZ() - pNode0->getZ();
    x = y1*z2 - z1*y2;
    y = z1*x2 - x1*z2;
    z = x1*y2 - y1*x2;
    return sqrt(x*x + y*y + z*z)*0.5;
}
double& CBoundaryFace::calcArea()
{
    Utility::CLogger* pLogger = Utility::CLogger::Instance();
    CNode *pNode0, *pNode1, *pNode2;
    switch(mnShapeType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        mArea = 0.0;
        pNode0 = mvBNode[0]->getNode();
        pNode1 = mvBNode[1]->getNode();
        pNode2 = mvBNode[2]->getNode();
        mArea += triArea(pNode0, pNode1, pNode2);
        pNode0 = mvBNode[2]->getNode();
        pNode1 = mvBNode[3]->getNode();
        pNode2 = mvBNode[0]->getNode();
        mArea += triArea(pNode0, pNode1, pNode2);
        break;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        pNode0 = mvBNode[0]->getNode();
        pNode1 = mvBNode[1]->getNode();
        pNode2 = mvBNode[2]->getNode();
        mArea = triArea(pNode0, pNode1, pNode2);
        break;
    default:
        pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::Area");
        break;
    }
    return mArea;
}
//--
// ディレクレ上位グリッドのEntValue分配 => 数式処理
//--
void CBoundaryFace::distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();
    uiint *pnEdgeVert;
    uiint iedge, numOfEdge;
    double dAveVal, calcVal;

    CHecMPI *pMPI= CHecMPI::Instance();//debug

    // 辺中心
    switch(mnShapeType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        numOfEdge=4;
        for(iedge=0; iedge < numOfEdge; iedge++) {
            pnEdgeVert = pEdgeTree->getQuadLocalNodeNum(iedge);
            dAveVal=0.0;
            ////    dAveVal += mvBNode[pnEdgeVert[0]]->getValue(dof, mgLevel);
            ////    dAveVal += mvBNode[pnEdgeVert[1]]->getValue(dof, mgLevel);
            dAveVal += mvBNode[pnEdgeVert[0]]->getEntValue(dof, mgLevel);
            dAveVal += mvBNode[pnEdgeVert[1]]->getEntValue(dof, mgLevel);
            dAveVal *= 0.5;
            if(mnShapeType==ElementType::Quad) {
                if(mgLevel!=nMaxMGLevel) {
                    ////    mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);
                    ////mvEdgeBNode[iedge]->initValue(dof, mgLevel+1);
                    mvEdgeBNode[iedge]->setEntValue(dof, mgLevel+1, dAveVal);//基礎データ

                    //////debug
                    ////cout << "BoundaryFace::distDirichletVal --- A  rank:" << pMPI->getRank() << endl;

                    calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
                    mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, calcVal);//-----Level+1: 数式処理あり=>数式計算値，数式なし=>基礎データ

                    //////debug
                    ////cout << "BoundaryFace::distDirichletVal --- B  dof:" << dof
                    ////        << "  rank:" << pMPI->getRank()
                    ////        << "  entVal:" << dAveVal
                    ////        << "  calcVal:" << calcVal << endl;
                }
            }
            if(mnShapeType==ElementType::Quad2) {
                ////    mvEdgeBNode[iedge]->setValue(dof, mgLevel,   dAveVal);
                ////mvEdgeBNode[iedge]->initValue(dof, mgLevel+1);
                mvEdgeBNode[iedge]->setEntValue(dof, mgLevel,   dAveVal);

                calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
                mvEdgeBNode[iedge]->setValue(dof, mgLevel, calcVal);//----- Level

                if(mgLevel!=nMaxMGLevel) {
                    ////   mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);
                    mvEdgeBNode[iedge]->setEntValue(dof, mgLevel+1, dAveVal);

                    calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
                    mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, calcVal);//----- Level+1
                }
            }
        };
        break;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        numOfEdge=3;
        for(iedge=0; iedge < numOfEdge; iedge++) {
            pnEdgeVert = pEdgeTree->getTriangleLocalNodeNum(iedge);
            dAveVal=0.0;
            ////    dAveVal += mvBNode[pnEdgeVert[0]]->getValue(dof, mgLevel);
            ////    dAveVal += mvBNode[pnEdgeVert[1]]->getValue(dof, mgLevel);
            dAveVal += mvBNode[pnEdgeVert[0]]->getEntValue(dof, mgLevel);
            dAveVal += mvBNode[pnEdgeVert[1]]->getEntValue(dof, mgLevel);
            dAveVal *= 0.5;
            if(mnShapeType==ElementType::Triangle) {
                if(mgLevel!=nMaxMGLevel) {
                    ////mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);
                    ////mvEdgeBNode[iedge]->initValue(dof, mgLevel+1);
                    mvEdgeBNode[iedge]->setEntValue(dof, mgLevel+1, dAveVal);

                    calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
                    mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, calcVal);//-----Level+1: 数式処理あり=>数式計算値，数式なし=>基礎データ
                }
            }
            if(mnShapeType==ElementType::Triangle2) {
                ////mvEdgeBNode[iedge]->setValue(dof, mgLevel,   dAveVal);
                ////mvEdgeBNode[iedge]->initValue(dof, mgLevel+1);
                mvEdgeBNode[iedge]->setEntValue(dof, mgLevel,   dAveVal);

                calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
                mvEdgeBNode[iedge]->setValue(dof, mgLevel, calcVal);//----- Level

                if(mgLevel!=nMaxMGLevel) {
                    ////mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);
                    mvEdgeBNode[iedge]->setEntValue(dof, mgLevel+1, dAveVal);

                    calcVal= getCalcValue( mvEdgeBNode[iedge], dAveVal, dof, vPolandDOF, mPoland );
                    mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, calcVal);//-----Level+1: 数式処理あり=>数式計算値，数式なし=>基礎データ
                }
            }
        };
        break;
    default:
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryFace::distDirichletVal, invalid ElementType");
        break;

    }//switch END

    uiint ivert, numOfVert=getNumOfVert();
    double dVal;
    dAveVal=0.0;
    // 頂点
    for(ivert=0; ivert < numOfVert; ivert++) {
        ////dVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dVal = mvBNode[ivert]->getEntValue(dof, mgLevel);
        dAveVal += dVal;
        if(mgLevel!=nMaxMGLevel) {
            ////mvBNode[ivert]->setValue(dof, mgLevel+1, dVal);
            ////mvBNode[ivert]->initValue(dof, mgLevel+1);
            mvBNode[ivert]->setEntValue(dof, mgLevel+1, dVal);

            calcVal= getCalcValue( mvBNode[ivert], dVal, dof, vPolandDOF, mPoland );
            mvBNode[ivert]->setValue(dof, mgLevel+1, calcVal);//-----Level+1: 数式処理あり=>数式計算値，数式なし=>基礎データ
        }
    };
    // 面中心
    dAveVal /= (double)numOfVert;
    if(mgLevel!=nMaxMGLevel) {
        ////mpFaceBNode->setValue(dof, mgLevel+1, dAveVal);
        ////mpFaceBNode->initValue(dof, mgLevel+1);
        mpFaceBNode->setEntValue(dof, mgLevel+1, dAveVal);

        calcVal= getCalcValue( mpFaceBNode, dAveVal, dof, vPolandDOF, mPoland );
        mpFaceBNode->setValue(dof, mgLevel+1, calcVal);//-----Level+1: 数式処理あり=>数式計算値，数式なし=>基礎データ
    }
}
void CBoundaryFace::replaceEdgeBNode()
{
    if(mnShapeType==ElementType::Quad2) {
        mvBNode.resize(NumberOfNode::Quad2());
        uiint nNumOfVert=NumberOfVertex::Quad();
        uiint nNumOfEdge=NumberOfEdge::Quad();
        uiint iedge;
        for(iedge=0; iedge < nNumOfEdge; iedge++) {
            CBoundaryNode *pBNode = mvEdgeBNode[iedge];
            mvBNode[nNumOfVert + iedge] = pBNode;
        };
    }
    if(mnShapeType==ElementType::Triangle2) {
        mvBNode.resize(NumberOfNode::Triangle2());
        uiint nNumOfVert=NumberOfVertex::Triangle();
        uiint nNumOfEdge=NumberOfEdge::Triangle();
        uiint iedge;
        for(iedge=0; iedge < nNumOfEdge; iedge++) {
            CBoundaryNode *pBNode = mvEdgeBNode[iedge];
            mvBNode[nNumOfVert + iedge] = pBNode;
        };
    }
}
void CBoundaryFace::deleteProgData()
{
    if(mnOrder==ElementOrder::First || mnOrder==ElementOrder::Zero) {
        vector<CBoundaryNode*>().swap(mvEdgeBNode);
    }
}
