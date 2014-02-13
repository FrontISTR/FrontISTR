/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryNodeMesh.h
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
#include "TypeDef.h"
#include "Logger.h"
#include "Node.h"
#include "BoundarySBNode.h"
#include <map>
#include "BoundaryType.h"

#include "Poland.h"

namespace pmw
{

#ifndef _BOUNDARYNODE_MESH_H
#define	_BOUNDARYNODE_MESH_H

class CBoundaryNodeMesh
{
public:
    CBoundaryNodeMesh();
    virtual ~CBoundaryNodeMesh();
protected:
    uiint mnID;
    uiint mnBndType;
    string msName;
    vector<CBoundarySBNode*>  mvBNode;
    map<uiint, uiint, less<uiint> > mmID2Index;
    map<uiint, uiint, less<uiint> > mmNodeID2BNodeID;

    vuint mvPolandDOF;//数式が定義されているDOF
    map<uiint, CPoland*> mmPoland;//自由度別 Poland

public:
    void setID(const uiint& boundID) {
        mnID= boundID;
    }
    uiint& getID() {
        return mnID;
    }
    void setName(const string& name) {
        msName = name;
    }
    string& getName() {
        return msName;
    }
    uiint getNumOfBNode() {
        return mvBNode.size();
    }
    void resizeBNode(const uiint& res_size);
    void setBNode(const uiint& index, CBoundarySBNode *pBNode);
    void addBNode(CBoundarySBNode *pBNode);
    CBoundarySBNode* getBNodeIX(const uiint& index) {
        return mvBNode[index];
    }
    CBoundarySBNode* getBNodeID(const uiint& id) {
        uiint index= mmID2Index[id];
        return mvBNode[index];
    }
    void  setBndType(const uiint& boundType);
    uiint& getBndType() {
        return mnBndType;
    }

    //--
    // Poland
    //--
    void setPoland(CPoland* pPoland, string sNumForm, const uiint& dof);//境界条件の数式セット
    ////void setBndFormula(string str, const uiint& dof);
    ////void execNumForm();//一括処理

    ////bool existPoland();
    bool existPoland(const uiint& dof);
    CPoland* getPoland(const uiint& dof);
};
#endif	/* _BOUNDARYNODE_MESH_H */
}
