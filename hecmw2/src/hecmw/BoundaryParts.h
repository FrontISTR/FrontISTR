/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BoundaryParts.h
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
    uiint mnElementID; 
    vector<CBoundaryNode*> mvBNode;
    map<uiint, uiint, less<uiint> > mmBNodeID2Index;
    uiint mnShapeType;
    uiint mnOrder;    
    map<uiint, double, less<uiint> > mmValue;
    CElement *mpElement;
public:
    void resizeBNode(const uiint& res_size);
    void setBNode(const uiint& ivert, CBoundaryNode *pBNode);
    uiint getNumOfBNode(){ return mvBNode.size();}
    CBoundaryNode* getBNode(const uiint& ibnode){ return mvBNode[ibnode];}
    uiint& getVertIndex(CBoundaryNode *pBNode){ return mmBNodeID2Index[pBNode->getID()];}
    uiint& getOrder(){ return mnOrder;}
    virtual uiint getNumOfVert()=0;
    void setID(const uiint& id){ mnID= id;}
    uiint& getID(){ return mnID;}
    void setElementID(const uiint& id){ mnElementID = id;}
    uiint& getElementID(){ return mnElementID;}
    void setElement(CElement *pElem){ mpElement= pElem;}
    CElement* getElement(){ return mpElement;}
    void setBndValue(const uiint& dof, const double& val){ mmValue[dof]= val;}
    double& getBndValue(const uiint& dof){ return mmValue[dof];}
    void setupVertexElemID();
};
#endif	/* _BOUNDARY_PARTS_H_ */
}
