/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   MeshFactory.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef MESH_FACTORY_HH_E74CAEBF_35A3_407b_A7B8_49D7A2DDDF8E
#define MESH_FACTORY_HH_E74CAEBF_35A3_407b_A7B8_49D7A2DDDF8E
#include "CommonStd.h"
#include "TypeDef.h"
#include "Logger.h"
#include "Element4Factory.h"
#include "ElementProperty.h"
#include "ScalarNode.h"
#include "VectorNode.h"
#include "ScalarVectorNode.h"
#include "NodeType.h"
#include "BoundaryNode.h"
#include "BoundaryFace.h"
#include "BoundaryVolume.h"
#include "BoundaryType.h"
#include "AggregateElement.h"
#include "AggregateNode.h"
#include "Mesh.h"
#include "AssyModel.h"
#include "GMGModel.h"
#include "IndexBucket.h"
#include "Material.h"
#include "CommunicationMesh.h"
#include "ProgElementTree.h"
#include "boost/lexical_cast.hpp"
#include "SkinFace.h"
#include "MasterFace.h"
#include "ContactMesh.h"
#include "ContactNode.h"
#include "CommMesh2.h"
#include "CommFace.h"
#include "CommNode.h"
namespace pmw{
class CMeshFactory
{
public:
    static CMeshFactory* Instance(){
        static CMeshFactory meshFactory;
        return &meshFactory;
    }
private:
    CMeshFactory(void);
public:
    virtual ~CMeshFactory(void);
protected:
    CGMGModel    *mpGMGModel;
    CAssyModel   *mpTAssyModel;
    CMesh        *mpTMesh;     
    CCommMesh    *mpTCommMesh; 
    CCommMesh2   *mpTCommMesh2;
    uint    mMGLevel;
    Utility::CLogger *mpLogger;
public:
    virtual void setGMGModel(CGMGModel* pGMGModel){ mpGMGModel = pGMGModel;}
    void GeneAssyModel(const uint& num_of_mgLevel);
    void reserveMesh(const uint& mgLevel,const uint& num_of_mesh);
    void GeneMesh(const uint& mgLevel, const uint& mesh_id, const uint& index);
    void setMGLevel(const uint& num_of_mglevel){ mMGLevel= num_of_mglevel;}
    uint& getMGLevel(){ return mMGLevel;}
    void refineMesh();
protected:
    void GeneProgElem(const uint& ilevel,CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
	void dividHexa(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividTetra(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividPrism(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividPyramid(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividQuad(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividTriangle(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividBeam(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
public:
    void reserveNode(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node);
    void GeneNode(const uint& mgLevel, const uint& mesh_id, const uint& id, const vdouble& coord,
                  const uint& nodeType, const uint& numOfScaParam, const uint& numOfVecParam);
    void setupNode(const uint& mgLevel, const uint& mesh_id);
    void reserveElement(const uint& mgLevel, const uint& mesh_id, const uint& num_of_element);
    void GeneElement(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& type, const vint& node_id);
    void setupElement(const uint& mgLevel, const uint& mesh_id);
    void reserveAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node);
    void GeneAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node);
    void reserveBoundaryNode(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void reserveBoundaryFace(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void reserveBoundaryVolume(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void GeneBoundaryNode(const uint& mgLevel, const uint& mesh_id, const uint& id,
                          const uint& dof, const uint& bndType,const vdouble& vVal);
    void GeneBoundaryFace(const uint& mgLevel, const uint& mesh_id, const uint& elem_id, const uint& face_id,
                          const uint& dof, const uint& bndType,const vdouble& vVal);
    void GeneBoundaryVolume(const uint& mgLevel, const uint& mesh_id, const uint& id,
                          const uint& dof, const uint& bndType,const vdouble& vVal);
    void initBucketNode(const uint& mgLevel, const uint& mesh_id, const uint& maxID, const uint& minID);
    void setIDBucketNode(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& index);
    void setupBucketNode(const uint& mgLevel, const uint& mesh_id);
    void initBucketElement(const uint& mgLevel, const uint& mesh_id, const uint& maxID, const uint& minID);
    void setIDBucketElement(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& index);
    void setupBucketElement(const uint& mgLevel, const uint& mesh_id);
    void setupBucketMesh(const uint& mgLevel, const uint& maxID, const uint& minID);
    void reserveMaterial(const uint& res_size);
    void GeneMaterial(const uint& mesh_id, const uint& material_id, string& name, vuint& vType, vdouble& vValue);
    void reserveCommMesh(const uint& mgLevel, const uint& mesh_id, const uint& res_size);
    void GeneCommMesh(const uint& mgLevel, const uint& mesh_id, const uint& comID, const uint& myRank, const uint& nTransmitRank);
    void reserveCommNode(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, const uint& res_size);
    void GeneCommNode(const uint& mgLevel, const uint& commNodeID, 
                      const uint& mesh_id, const uint& commesh_id, const uint& nodeID, const uint& rank);
    void reserveCommElement(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, const uint& res_size);
    void GeneCommElement(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, 
                                   const uint& nType, const uint& elemID, vuint& vCommNodeID);
protected:
    void GeneProgCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem);
    void dividCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem);
public:
    void GeneContactMesh(const uint& contactID, const uint& myRank, const uint& tranRank);
    void GeneContactNode(const uint& mgLevel, const uint& contactID, const uint& conNodeID, const vdouble& vCoord,
            const string& s_param_type, const uint& numOfVector, const uint& numOfScalar,
            bool bmesh, const uint& meshID, const uint& nodeID,
            const uint& rank, const uint& maslave);
    void GeneMasterFace(const uint& contactID, const uint& shapeType, const uint& masterFaceID,
            bool bmesh, const uint& meshID, const uint& elemID, const uint& elemFaceID,
            const vuint& vConNodeID, const uint& face_rank);
    void GeneSlaveFace(const uint& contactID, const uint& shapeType, const uint& slaveFaceID,
            bool bmesh, const uint& meshID, const uint& elemID, const uint& elemFaceID,
            const vuint& vConNodeID, const uint& face_rank);
    void refineContactMesh();
protected:
    void setProgHexaMPCMaster(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l);
    void setProgHexaMPCSlave(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l);
    void setProgHexaCommMesh2Entity(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l);
public:
    void refineCommMesh2();
    void GeneCommMesh2(const uint& mgLevel, const uint& mesh_id, const uint& comID, const uint& numOfFace, const uint& numOfCommNode,
            const uint& myRank, const uint& nTransmitRank);
    void GeneCommFace(const uint& mgLevel, const uint& commeshID, const uint& face_id,
            const uint& mesh_id,const uint elem_id, const uint& elem_ent_num, const vuint& vCommNodeID);
    void GeneCommNodeCM2(const uint& mgLevel, const uint& mesh_id, const uint& node_id, const uint& commeshID, const uint& comm_node_id, const vdouble& vCoord);
};
}
#endif
