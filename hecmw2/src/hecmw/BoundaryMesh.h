/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BoundaryMesh.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "TypeDef.h"
#include "BoundaryNode.h"
#include <map>
namespace pmw{
#ifndef _BOUNDARYMESH_H
#define	_BOUNDARYMESH_H
class CBoundaryMesh{
public:
    CBoundaryMesh();
    virtual ~CBoundaryMesh();
protected:
    uiint mnID;     
    uiint mnBndType;
    string msName; 
    uiint mMGLevel;   
    uiint mMaxMGLevel;
    vector<CBoundaryNode*> mvBNode;
    map<uiint, uiint, less<uiint> > mmBNodeID2Index;
    vuint  mvDOF;                            
    map<uiint, uiint, less<uiint> > mmDOF2Index;
    uiint mnEdgeNodeCount;
public:
    void setID(const uiint& id){ mnID= id;}
    uiint& getID(){ return mnID;}
    void setName(const string& name){ msName = name;}
    string& getName(){ return msName;}
    void addDOF(const uiint& dof);
    void setDOF(const uiint& index, const uiint& dof);
    void resizeDOF(const uiint& res_size);
    uiint& getDOF(const uiint& index);
    uiint& getDOF_Index(const uiint& dof);
    uiint getNumOfDOF();
    void setBndType(const uiint& bndType){ mnBndType= bndType;}
    uiint& getBndType(){ return mnBndType;}
    void setMGLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getMGLevel(){ return mMGLevel;}
    void setMaxMGLevel(const uiint& maxLevel){ mMaxMGLevel= maxLevel;}
    uiint& getMaxMGLevel(){ return mMaxMGLevel;}
    void resizeCGrid_BNodeValue(const uiint& maxLevel);
    void resizeBNode(const uiint& res_size){ mvBNode.resize(res_size);}
    void setBNode(const uiint& index, CBoundaryNode *pBNode);
    void addBNode(CBoundaryNode* pBNode);
    uiint getNumOfBNode(){ return mvBNode.size();}
    CBoundaryNode* getBNodeIX(const uiint& index){ return mvBNode[index];}
    CBoundaryNode* getBNodeID(const uiint& id){
        uiint index= mmBNodeID2Index[id];
        return mvBNode[index];
    }
    uiint& getBNodeIndex(const uiint& id){ return mmBNodeID2Index[id];}
protected:
    virtual void distNeumannValue()=0;
    virtual void distDirichletValue()=0;
public:
    void distValueBNode();
};
#endif	/* _BOUNDARYMESH_H */
}
