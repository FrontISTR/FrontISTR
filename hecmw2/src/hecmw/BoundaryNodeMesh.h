/*
 *
 * File:   BoundaryNodeMesh.h
 * Author: ktakeda
 *
 * BoundaryNodeを面として管理(MultiGrid対応)
 *
 *
 * Created on 2010/04/07, 14:37
 */
#include "TypeDef.h"
#include "Logger.h"

#include "Node.h"
//#include "BoundaryNode.h"
#include "BoundarySBNode.h"//BoundaryNodeMesh用 BNode
#include <map>// ID->index

#include "BoundaryType.h"

namespace pmw{
#ifndef _BOUNDARYNODE_MESH_H
#define	_BOUNDARYNODE_MESH_H
class CBoundaryNodeMesh{
public:
    CBoundaryNodeMesh();
    virtual ~CBoundaryNodeMesh();

protected:
    uint mnID;  //BoundaryID
    uint mnBndType;//Dirichlet | Neumann

    string msName; //境界名称

    // *節点への境界条件は,全ての階層で同じ*
    // *全てのMGLevelで同一のBoundaryNodeMeshを使用*

    //  mnDOF 境界条件の自由度Index番号 => BNodeへ : 節点の境界自由度は節点毎にバラバラで対応

    vector<CBoundarySBNode*>  mvBNode;//境界節点

    map<uint, uint, less<uint> > mmID2Index; // ID -> Index
    map<uint, uint, less<uint> > mmNodeID2BNodeID; // NodeID -> BoundaryNodeID


public:
    // BoundaryID
    void setID(const uint& boundID){ mnID= boundID;}
    uint& getID(){ return mnID;}

    //境界名称
    void setName(const string& name){ msName = name;}
    string& getName(){ return msName;}


    uint getNumOfBNode(){ return mvBNode.size();}
    void resizeBNode(const uint& res_size);
    void setBNode(const uint& index, CBoundarySBNode *pBNode);
    void addBNode(CBoundarySBNode *pBNode);
    CBoundarySBNode* getBNodeIX(const uint& index){ return mvBNode[index];}
    CBoundarySBNode* getBNodeID(const uint& id){
        uint index= mmID2Index[id];
        return mvBNode[index];
    }


    // Dirichlet | Neumann
    void  setBndType(const uint& boundType);
    uint& getBndType(){ return mnBndType;}

    
};
#endif	/* _BOUNDARYNODE_MESH_H */
}



