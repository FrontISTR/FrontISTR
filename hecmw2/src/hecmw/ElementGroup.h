/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ElementGroup.h
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
#include "CommonStd.h"
#include "TypeDef.h"
namespace pmw
{
class CMesh;
class CElement;
#ifndef ELEMENTGROUP_H
#define	ELEMENTGROUP_H
class CElementGroup
{
public:
    CElementGroup();
    virtual ~CElementGroup();
protected:
    CMesh *mpMesh;
    vuint mvElementID;
    uiint mnGrpID;
    string msGrpName;
public:
    void setMesh(CMesh *pMesh);
    void setID(const uiint& id);
    uiint& getID();
    void setName(const string& sGrpName);
    string& getName();
    uiint getNameLength();
    void addElementID(const uiint& nElemID);
    uiint getNumOfElementID();
    uiint& getElementID(const uiint& index);
    void refine(CElementGroup* pProgElemG);
};
#endif	/* ELEMENTGROUP_H */
}
