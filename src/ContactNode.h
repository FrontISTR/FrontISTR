/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ContactNode.h
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
#include "Vertex.h"
#include <vector>
#include <iostream>
namespace pmw{
#ifndef _CONTACTNODE_H
#define	_CONTACTNODE_H
class CContactNode:public CVertex{
public:
    CContactNode();
    virtual ~CContactNode();
protected:
    uint myRank; 
    uint mMeshID;
    uint mNodeID;
    uint mLevel;
    bool mbMesh;
    bool mbNode;
    vdouble mvDisplacement;
    vdouble mvScalar;
    uint mMFaceID;
    bool mbSlave; 
    bool mbMarkingMFace;
    map<uint, uint, less<uint> > mmMasterFaceID;
    vector<bool>                mvbMasterFaceID;
    vector<uint> mvKnotID;
public:
    void markingSelfMesh();
    void markingSelfNode();
	void  setLevel(const uint& level){ mLevel= level;}
    uint& getLevel(){ return mLevel;}
	void  pushLevelMarking(){ mvbMasterFaceID.push_back(false);}
    void  setRank(const uint& rank){ myRank= rank;}
    uint& getRank(){ return myRank;}
    void  setMeshID(const uint& meshID){ mMeshID= meshID;}
    uint& getMeshID(){ return mMeshID;}
    void  setNodeID(const uint& nodeID){ mNodeID= nodeID;}
    uint& getNodeID(){ return mNodeID;}
    void resizeDisp(const uint& dof);
    void initDisp();
    void setDisp(const uint& idof, const double& disp){ mvDisplacement[idof]= disp;}
    double& getDisp(const uint& idof){ return mvDisplacement[idof];}
    uint getNumOfDisp(){return mvDisplacement.size();}
    void resizeScalar(const uint& numOfScalar);
    void initScalar();
    void setScalar(const uint& i, const double& val){ mvScalar[i]=val;}
    double& getScalar(const uint& i){ return mvScalar[i];}
    uint getNumOfScalar(){return mvScalar.size();}
    void markingSlave();
    bool isSlave(){ return mbSlave;}
    void setMasterFaceID(const uint& faceID, const uint& level);
    uint& getMasterFaceID(const uint& level);
    bool have_MasterFaceID(const uint& level);
    void  resizeOctreeID(const uint& res_size);
    void  setOctreeID(const uint& layer, const uint& knot_id);
    uint& getOctreeID(const uint& layer){ return mvKnotID[layer];}
};
#endif	/* _CONTACTNODE_H */
}
