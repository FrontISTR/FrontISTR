/* 
 * File:   BoundaryParts.h
 * Author: ktakeda
 *
 * Created on 2010/06/09, 14:59
 */
#include "TypeDef.h"
#include <map>

#include "BoundaryNode.h"
#include "Logger.h"


#include "Node.h"
#include "Element.h"

namespace pmw{
#ifndef _BOUNDARY_PARTS_H_
#define	_BOUNDARY_PARTS_H_
class CBoundaryParts{
public:
    CBoundaryParts();
    virtual ~CBoundaryParts();

protected:
    uint mnID;
    uint mnElementID; //Mesh-Element ID

    vector<CBoundaryNode*> mvBNode;// 頂点(値を分配するNode)
    map<uint, uint, less<uint> > mmBNodeID2Index;

    uint mnShapeType;// Beam, Quad, Triangle, Tetra, Hexa
    uint mnOrder;    // 1次 .or. 2次
    
    ////vuint mvDOF;// => BounaryMeshに移行.
    map<uint, double, less<uint> > mmValue;//Partsが所有するDOF境界値


    CElement *mpElement;//Boundaryが載ってる要素

public:
    // 頂点(境界Node)
    void resizeBNode(const uint& res_size);
    void setBNode(const uint& ivert, CBoundaryNode *pBNode);
    uint getNumOfBNode(){ return mvBNode.size();}
    CBoundaryNode* getBNode(const uint& ivert){ return mvBNode[ivert];}

    uint& getVertIndex(CBoundaryNode *pBNode){ return mmBNodeID2Index[pBNode->getID()];}

    // 1次 .or. 2次
    uint& getOrder(){ return mnOrder;}

    // 頂点数
    virtual uint getNumOfVert()=0;

    // 境界パーツID(BFaceID,BEdgeID,BVolumeID)
    void setID(const uint& id){ mnID= id;}
    uint& getID(){ return mnID;}

    // Mesh-Element-ID(Element ID)
    void setElementID(const uint& id){ mnElementID = id;}
    uint& getElementID(){ return mnElementID;}

    // Mesh-Element
    void setElement(CElement *pElem){ mpElement= pElem;}
    CElement* getElement(){ return mpElement;}

    // 境界値
    void setBndValue(const uint& dof, const double& val){ mmValue[dof]= val;}
    double& getBndValue(const uint& dof){ return mmValue[dof];}

    
    // ----
    // 頂点に(Edge,Face,Volume)IDの集合をセット
    // ----
    void setupVertexElemID();
};
#endif	/* _BOUNDARY_PARTS_H_ */
}

