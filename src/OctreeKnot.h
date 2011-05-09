/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   OctreeKnot.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "TypeDef.h"
#include "ContactNode.h"
#include <iostream>
namespace pmw{
#ifndef _OCTREEKNOT_H
#define	_OCTREEKNOT_H
class COctreeKnot{
public:
    COctreeKnot();
    virtual ~COctreeKnot();
protected:
    uint mLayer;
    uint mID;   
    COctreeKnot* mpParentKnot;       
    vector<COctreeKnot*> mvChildKnot;
    vector<CContactNode*> mvMasterNode;
    vector<CContactNode*> mvSlaveNode; 
    double minX,maxX;
    double minY,maxY;
    double minZ,maxZ;
public:
    void setLayerID(const uint& layer){ mLayer= layer;}
    void setID(const uint& id){ mID= id;}
    uint& getLayerID(){ return mLayer;}
    uint& getID(){ return mID;}
    void setItemProp();
    void reserveMasterNode(const uint& res_size){ mvMasterNode.reserve(res_size);}
    void addMasterNode(CContactNode* pConNode){ mvMasterNode.push_back(pConNode);}
    uint getNumOfMasterNode(){ return mvMasterNode.size();}
    vector<CContactNode*>& getMasterNode(){ return mvMasterNode;}
    CContactNode*  getMasterNode(const uint& index){ return mvMasterNode[index];}
    void reserveSlaveNode(const uint& res_size){ mvSlaveNode.reserve(res_size);}
    void addSlaveNode(CContactNode* pConNode){ mvSlaveNode.push_back(pConNode);}
    uint getNumOfSlaveNode(){ return mvSlaveNode.size();}
    vector<CContactNode*>& getSlaveNode(){ return mvSlaveNode;}
    CContactNode*  getSlaveNode(const uint& index){ return mvSlaveNode[index];}
    void setX(const double& min, const double& max){ minX=min; maxX=max;}
    void setY(const double& min, const double& max){ minY=min; maxY=max;}
    void setZ(const double& min, const double& max){ minZ=min; maxZ=max;}
    double& getMinX(){return minX;} double& getMaxX(){return maxX;}
    double& getMinY(){return minY;} double& getMaxY(){return maxY;}
    double& getMinZ(){return minZ;} double& getMaxZ(){return maxZ;}
    void setParentKnot(COctreeKnot* pKnot){ mpParentKnot= pKnot;}
    COctreeKnot* getParentKnot(){ return mpParentKnot;}
    void createChildKnot();
    COctreeKnot* getChildKnot(const uint& pos){ return mvChildKnot[pos];}
    void distItem();
};
#endif	/* _OCTREEKNOT_H */
}
