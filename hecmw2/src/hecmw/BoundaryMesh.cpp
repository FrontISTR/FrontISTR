/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryMesh.cpp
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
#include "BoundaryMesh.h"
using namespace pmw;
CBoundaryMesh::CBoundaryMesh()
{
    mnEdgeNodeCount = 0;
}
CBoundaryMesh::~CBoundaryMesh()
{
    if(mMaxMGLevel==mMGLevel) for_each(mvBNode.begin(), mvBNode.end(), DeleteObject());

    if(mMaxMGLevel==mMGLevel) {
        map<uiint,CPoland*>::iterator it;
        for(it=mmPoland.begin(); it != mmPoland.end(); it++) {
            CPoland* pPoland= it->second;
            delete pPoland;
        };
    }
}
void CBoundaryMesh::addDOF(const uiint& dof)
{
    mvDOF.push_back(dof);
    mmDOF2Index[dof] = mvDOF.size() - 1;
}
void CBoundaryMesh::setDOF(const uiint& index, const uiint& dof)
{
    mvDOF[index] = dof;
    mmDOF2Index[dof] = index;
}
void CBoundaryMesh::resizeDOF(const uiint& res_size)
{
    mvDOF.resize(res_size);
}
uiint& CBoundaryMesh::getDOF(const uiint& index)
{
    return mvDOF[index];
}
uiint& CBoundaryMesh::getDOF_Index(const uiint& dof)
{
    return mmDOF2Index[dof];
}
uiint CBoundaryMesh::getNumOfDOF()
{
    return mvDOF.size();
}
//--
// 制御ファイルでリファイン数を指定するので、リファイン関数でコースグリッドのBNodeをresize
// # メッシュ・ファイル入力時には、Level数が不明なので1個だけ確保しておき、後からresizeをする.
//--
void CBoundaryMesh::resizeCGrid_BNodeValue(const uiint& maxLevel)
{
    uiint i, nNumOfBNode= mvBNode.size();
    for(i=0; i < nNumOfBNode; i++) {
        mvBNode[i]->resizeValue(maxLevel+1);
    };
}
void CBoundaryMesh::setBNode(const uiint& index, CBoundaryNode* pBNode)
{
    uiint id;
    id= pBNode->getID();
    mmBNodeID2Index[id]= index;
    mvBNode[index]= pBNode;
}
void CBoundaryMesh::addBNode(CBoundaryNode* pBNode)
{
    mvBNode.push_back(pBNode);
    uiint id;
    id= pBNode->getID();
    mmBNodeID2Index[id]= mvBNode.size()-1;
}
void CBoundaryMesh::distValueBNode()
{
    switch(mnBndType) {
    case(BoundaryType::Neumann):
        distNeumannValue();
        break;
    case(BoundaryType::Dirichlet):
        distDirichletValue();
        break;
    }
}

//--
// Poland
//--
void CBoundaryMesh::setPoland(CPoland* pPoland, string sNumForm, const uiint& dof)
{
    mvPolandDOF.push_back(dof);
    stable_sort( mvPolandDOF.begin(), mvPolandDOF.end() );

    mmPoland[dof]=pPoland;

    mmPoland[dof]->setOrigForm(sNumForm);
    mmPoland[dof]->transForm();
}

CPoland* CBoundaryMesh::getPoland(const uiint& dof)
{
    return mmPoland[dof];
}
map<uiint,CPoland*>& CBoundaryMesh::getPoland()
{
    return mmPoland;
}
vuint& CBoundaryMesh::getPolandDOF()
{
    return mvPolandDOF;
}
////bool CBoundaryMesh::existPoland()
////{
////    if( mmPoland.empty() ){
////        return false;
////    }else{
////        return true;
////    }
////}
bool CBoundaryMesh::existPoland(const uiint& dof)
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

void CBoundaryMesh::setPoland(map<uiint,CPoland*> mPoland)
{
    mmPoland= mPoland;
}
void CBoundaryMesh::setPolandDOF(vuint& vPolandDOF)
{
    mvPolandDOF= vPolandDOF;
}


