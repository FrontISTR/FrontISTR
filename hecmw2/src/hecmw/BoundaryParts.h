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
    uiint mnID;
    uiint mnElementID; //Mesh-Element ID

    vector<CBoundaryNode*> mvBNode;// 頂点(値を分配するNode)
    map<uiint, uiint, less<uiint> > mmBNodeID2Index;

    uiint mnShapeType;// Beam, Quad, Triangle, Tetra, Hexa
    uiint mnOrder;    // 1次 .or. 2次
    
    ////vuint mvDOF;// => BounaryMeshに移行.
    map<uiint, double, less<uiint> > mmValue;//Partsが所有するDOF境界値


    CElement *mpElement;//Boundaryが載ってる要素

public:
    // 頂点(境界Node)
    void resizeBNode(const uiint& res_size);
    void setBNode(const uiint& ivert, CBoundaryNode *pBNode);
    uiint getNumOfBNode(){ return mvBNode.size();}
    CBoundaryNode* getBNode(const uiint& ibnode){ return mvBNode[ibnode];}

    uiint& getVertIndex(CBoundaryNode *pBNode){ return mmBNodeID2Index[pBNode->getID()];}

    // 1次 .or. 2次
    uiint& getOrder(){ return mnOrder;}

    // 頂点数
    virtual uiint getNumOfVert()=0;

    // 境界パーツID(BFaceID,BEdgeID,BVolumeID)
    void setID(const uiint& id){ mnID= id;}
    uiint& getID(){ return mnID;}

    // Mesh-Element-ID(Element ID)
    void setElementID(const uiint& id){ mnElementID = id;}
    uiint& getElementID(){ return mnElementID;}

    // Mesh-Element
    void setElement(CElement *pElem){ mpElement= pElem;}
    CElement* getElement(){ return mpElement;}

    // 境界値
    void setBndValue(const uiint& dof, const double& val){ mmValue[dof]= val;}
    double& getBndValue(const uiint& dof){ return mmValue[dof];}

    
    // ----
    // 頂点に(Edge,Face,Volume)IDの集合をセット
    // ----
    void setupVertexElemID();
};
#endif	/* _BOUNDARY_PARTS_H_ */
}

