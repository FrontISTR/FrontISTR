//
// ElementGroup.cpp
//
//              2010.10.14
//              k.Takeda
#include "ElementGroup.h"
#include "Mesh.h"
using namespace pmw;


CElementGroup::CElementGroup()
{
    ;
}
CElementGroup::~CElementGroup()
{
    ;
}


// ElemGrp自身が所属しているMeshポインター
//
void CElementGroup::setMesh(CMesh* pMesh)
{
    mpMesh = pMesh;
}

// ID
//
void CElementGroup::setID(const uint& id)
{
    mnGrpID = id;
}
uint& CElementGroup::getID()
{
    return mnGrpID;
}

// Name
//
void CElementGroup::setName(const string& sGrpName)
{
    msGrpName = sGrpName;
}
string& CElementGroup::getName()
{
    return msGrpName;
}
uint CElementGroup::getNameLength()
{
    return (uint)msGrpName.length();
}

// ElementGrpへElemIDを追加
//
void CElementGroup::addElementID(const uint& nElemID)
{
    mvElementID.push_back(nElemID);
}

// ElementGrpのElemID数
//
uint CElementGroup::getNumOfElementID()
{
    return mvElementID.size();
}

// ElementGrpのElemID
//
uint& CElementGroup::getElementID(const uint& index)
{
    return mvElementID[index];
}


// 上位GridのElementGroupを生成
//
void CElementGroup::refine(CElementGroup* pProgElemG)
{
    uint iElem, nNumOfElem = mvElementID.size();
    uint nElemID;
    CElement *pElem, *pProgElem;
    for(iElem=0; iElem < nNumOfElem; iElem++){
        nElemID = mvElementID[iElem];
        pElem = mpMesh->getElement(nElemID);

        //cout << "ElementGroup::refine, pElem ID = " << pElem->getID() << endl;

        uint ivert, nNumOfVert = pElem->getNumOfNode();
        for(ivert=0; ivert < nNumOfVert; ivert++){
            pProgElem = pElem->getProgElem(ivert);

            pProgElemG->addElementID(pProgElem->getID());//// 上位GridのElementGrpへElemIDを追加.
        };
    };
}







