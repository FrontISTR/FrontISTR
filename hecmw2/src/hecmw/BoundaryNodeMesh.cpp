/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryNodeMesh.cpp
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
#include "BoundaryNodeMesh.h"
#include "Calc.h"
using namespace pmw;
CBoundaryNodeMesh::CBoundaryNodeMesh()
{
    ;
}
CBoundaryNodeMesh::~CBoundaryNodeMesh()
{
    map<uiint,CPoland*>::iterator it;
    for(it=mmPoland.begin(); it != mmPoland.end(); it++) {
        CPoland *pPoland= it->second;
        delete pPoland;
    };
}
void CBoundaryNodeMesh::setBndType(const uiint& boundType)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    switch(boundType) {
    case(BoundaryType::Dirichlet):
        mnBndType= boundType;
        break;
    case(BoundaryType::Neumann):
        mnBndType= boundType;
        break;
    default:
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error, CBoundaryNodeMesh::setType");
        break;
    }
}
void CBoundaryNodeMesh::resizeBNode(const uiint& res_size)
{
    mvBNode.resize(res_size);
}
void CBoundaryNodeMesh::setBNode(const uiint& index, CBoundarySBNode* pBNode)
{
    uiint id= pBNode->getID();
    mmID2Index[id]= index;
    CNode *pNode= pBNode->getNode();
    uiint node_id= pNode->getID();
    mmNodeID2BNodeID[node_id]= pBNode->getID();
}
void CBoundaryNodeMesh::addBNode(CBoundarySBNode* pBNode)
{
    mvBNode.push_back(pBNode);
    uiint index= mvBNode.size()-1;
    uiint id= pBNode->getID();
    mmID2Index[id]= index;
    CNode *pNode= pBNode->getNode();
    uiint node_id= pNode->getID();
    mmNodeID2BNodeID[node_id]= pBNode->getID();
}


//--
// Poland
//--
void CBoundaryNodeMesh::setPoland(CPoland* pPoland, string sNumForm, const uiint& dof)
{
    mvPolandDOF.push_back(dof);
    stable_sort( mvPolandDOF.begin(),mvPolandDOF.end() );

    mmPoland[dof]=pPoland;

    mmPoland[dof]->setOrigForm(sNumForm);
    mmPoland[dof]->transForm();
}
////void CBoundaryNodeMesh::setBndFormula(string str, const uiint& dof)
////{
////    mmPoland[dof]->setOrigForm(str);
////    mmPoland[dof]->transForm();
////}
////
////void CBoundaryNodeMesh::execNumForm()
////{
////    CCalc *pCalc= CCalc::Instance();
////
////    CBoundarySBNode *pBNode;
////    for(uiint i=0; i < mvBNode.size(); i++){
////        pBNode=mvBNode[i];
////        for(uiint idof=0; idof < pBNode->getNumOfDOF(); idof++){
////            uiint dof= pBNode->getDOF(idof);
////            double entVal=pBNode->getEntValue(dof);
////            double x=pBNode->getX(), y=pBNode->getY(), z=pBNode->getZ();
////            double val;
////            if(!mmPoland.empty()){
////                //数式処理
////                pCalc->setElementParam(entVal, x, y, z);
////                val= pCalc->Exec(mmPoland[dof]);
////
////                pBNode->setValue(dof, val);
////            }else{
////                //数式が無いのでentValが境界値
////                pBNode->setValue(dof, entVal);
////            }
////        };
////    };
////}


////bool CBoundaryNodeMesh::existPoland()
////{
////    if(!mmPoland.empty()){
////        return true;
////    }else{
////        return false;
////    }
////}

bool CBoundaryNodeMesh::existPoland(const uiint& dof)
{
    if(!mvPolandDOF.empty()) {
        uiint nNumOfDOF=mvPolandDOF.size();
        uiint low = 0;
        uiint high= nNumOfDOF-1;

        while(low <= high) {
            uiint ix=(low+high)/2;
            if(dof==mvPolandDOF[ix]) {
                return true;
            } else if(dof < mvPolandDOF[ix]) {
                if(ix!=0) {
                    high= ix-1;
                } else {
                    return false;
                }
            } else {
                low= ix+1;
            }
        };
        return false;
    } else {
        return false;
    }
}

CPoland* CBoundaryNodeMesh::getPoland(const uiint& dof)
{
    return mmPoland[dof];
}





