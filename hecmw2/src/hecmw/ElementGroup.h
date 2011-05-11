/* 
 * File:   ElementGroup.h
 * Author: ktakeda
 *
 * Created on 2010/10/14, 17:19
 */
#include "CommonStd.h"
#include "TypeDef.h"

namespace pmw{

class CMesh;
class CElement;

#ifndef ELEMENTGROUP_H
#define	ELEMENTGROUP_H
class CElementGroup{
public:
    CElementGroup();
    virtual ~CElementGroup();

protected:
    CMesh *mpMesh;
    vuint mvElementID;

    uiint mnGrpID;
    string msGrpName;

public:
    void setMesh(CMesh *pMesh);//ElemGrp自身が所属しているMeshポインター

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





