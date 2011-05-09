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

    uint mnGrpID;
    string msGrpName;

public:
    void setMesh(CMesh *pMesh);//ElemGrp自身が所属しているMeshポインター

    void setID(const uint& id);
    uint& getID();

    void setName(const string& sGrpName);
    string& getName();
    uint getNameLength();

    void addElementID(const uint& nElemID);

    uint getNumOfElementID();
    uint& getElementID(const uint& index);

    void refine(CElementGroup* pProgElemG);
};
#endif	/* ELEMENTGROUP_H */
}





