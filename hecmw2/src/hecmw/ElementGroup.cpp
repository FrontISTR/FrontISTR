/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ElementGroup.cpp
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
void CElementGroup::setMesh(CMesh* pMesh)
{
    mpMesh = pMesh;
}
void CElementGroup::setID(const uiint& id)
{
    mnGrpID = id;
}
uiint& CElementGroup::getID()
{
    return mnGrpID;
}
void CElementGroup::setName(const string& sGrpName)
{
    msGrpName = sGrpName;
}
string& CElementGroup::getName()
{
    return msGrpName;
}
uiint CElementGroup::getNameLength()
{
    return (uiint)msGrpName.length();
}
void CElementGroup::addElementID(const uiint& nElemID)
{
    mvElementID.push_back(nElemID);
}
uiint CElementGroup::getNumOfElementID()
{
    return mvElementID.size();
}
uiint& CElementGroup::getElementID(const uiint& index)
{
    return mvElementID[index];
}
void CElementGroup::refine(CElementGroup* pProgElemG)
{
    uiint iElem, nNumOfElem = mvElementID.size();
    uiint nElemID;
    CElement *pElem, *pProgElem;
    for(iElem=0; iElem < nNumOfElem; iElem++) {

        nElemID = mvElementID[iElem];
        pElem = mpMesh->getElement(nElemID);

        uiint ivert, nNumOfVert = pElem->getNumOfVert();

        for(ivert=0; ivert < nNumOfVert; ivert++) {
            pProgElem = pElem->getProgElem(ivert);
            pProgElemG->addElementID(pProgElem->getID());
        };
    };
}
