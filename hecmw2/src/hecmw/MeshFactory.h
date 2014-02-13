/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/MeshFactory.h
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
#include "CommonStd.h"
#include "TypeDef.h"
#include "Logger.h"
#include "Element4Factory.h"
#include "ElementProperty.h"
#include "ScalarNode.h"
#include "VectorNode.h"
#include "ScalarVectorNode.h"
#include "NodeType.h"
#include "BoundaryNodeMesh.h"
#include "BoundaryFaceMesh.h"
#include "BoundaryVolumeMesh.h"
#include "BoundaryEdgeMesh.h"
#include "BoundaryNode.h"
#include "BoundarySBNode.h"
#include "BoundaryFace.h"
#include "BoundaryHexa.h"
#include "BoundaryTetra.h"
#include "BoundaryPrism.h"
#include "BoundaryEdge.h"
#include "BoundaryType.h"
#include "BNodeMeshGrp.h"
#include "AggregateElement.h"
#include "AggregateNode.h"
#include "Mesh.h"
#include "AssyModel.h"
#include "GMGModel.h"
#include "IndexBucket.h"
#include "Material.h"
#include "CommunicationMesh.h"
#include "ProgElementTree.h"

#include "SkinFace.h"
#include "MasterFace.h"
#include "ContactMesh.h"
#include "ContactNode.h"
#include "CommMesh2.h"
#include "CommFace.h"
#include "CommNode.h"
#include "SolutionType.h"
#include "ElementGroup.h"
#include "HEC_MPI.h"
namespace pmw
{
#ifndef _MESH_FACTORY_HH_E74C
#define _MESH_FACTORY_HH_E74C
class CMeshFactory
{
public:
    static CMeshFactory* Instance() {
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
    uiint    mMGLevel;
    uiint mnSolutionType;
    Utility::CLogger *mpLogger;
public:
    virtual void setGMGModel(CGMGModel* pGMGModel) {
        mpGMGModel = pGMGModel;
    }
    void GeneAssyModel(const uiint& nNumOfLevel);
    void setGlobalCommData(const uiint& mgLevel, const uiint& nNumGlobalComm, const vector<pair<uiint, uiint> >& vCommPair, const vuint& vMeshID_CommID);
    void reserveMesh(const uiint& mgLevel,const uiint& num_of_mesh);
    void GeneMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& index, const uiint& nProp);
    void setSolutionType(const uiint& nSolutionType) {
        mnSolutionType = nSolutionType;
    }
    void setMGLevel(const uiint& num_of_mglevel);
    uiint& getMGLevel() {
        return mMGLevel;
    }
    void setupNodeGridSize();

    void refineMesh();

protected:
    void MGMeshConstruct();
    void SGMeshConstruct();
    void GeneProgElem(const uiint& ilevel,CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividHexa(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividTetra(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividPrism(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividPyramid(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividQuad(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividTriangle(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividBeam(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
public:
    void reserveNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_node);
    void GeneNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const vdouble& coord,
                  const uiint& nodeType, const uiint& nNumOfSDOF, const uiint& nNumOfVDOF);
    void setupNode(const uiint& mgLevel, const uiint& mesh_id);
    void reserveElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_element);
    void GeneElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const uiint& type, const vint& node_id);
    void setupElement(const uiint& mgLevel, const uiint& mesh_id);
    void resizeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_node);
    void GeneAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_node);

    void reserveBoundaryNodeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd);
    void GeneBoundaryNodeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name, map<uiint,string> mStrForm);
    void reserveBoundaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd);
    void GeneBoundaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name,
                              const uiint& numOfDOF, const vuint& vDOF, map<uiint,string> mStrForm);
    void reserveBoundaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd);
    void GeneBoundaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name,
                                const uiint& numOfDOF, const vuint& vDOF, map<uiint,string> mStrForm);
    void reserveBoundaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd);
    void GeneBoundaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name,
                              const uiint& numOfDOF, const vuint& vDOF, map<uiint,string> mStrForm);
    uiint getNumOfBounaryNodeMesh(const uiint& mgLevel, const uiint& mesh_id);
    uiint getNumOfBounaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id);
    uiint getNumOfBounaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id);
    uiint getNumOfBounaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id);
    CBoundaryFaceMesh* getBoundaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    CBoundaryVolumeMesh* getBoundaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    CBoundaryEdgeMesh* getBoundaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void GeneBoundaryNode(const uiint& mgLevel, const uiint& bnd_id,  const uiint& bndType,
                          const uiint& mesh_id, const uiint& node_id,
                          const uiint& b_node_id, const uiint& dof, const double& val);
    void initFaceAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void resizeFaceAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void GeneBoundaryFaceNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
                              const uiint& mesh_id, const uiint& node_id,
                              const uiint& b_node_id);
    void setValue_BoundaryFaceNode(const uiint& mesh_id, const uiint& bnd_id, const uiint& bnode_id, const uiint& dof, const double& val);
    void GeneBoundaryFace(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType, const uiint& elemType,
                          const uiint& mesh_id, const uiint& elem_id, const uiint& face_id, vuint& vBNodeID,
                          const uiint& b_face_id, const uiint& dof, const double& val);
    void initVolumeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void resizeVolumeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void GeneBoundaryVolumeNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
                                const uiint& mesh_id, const uiint& node_id,
                                const uiint& b_node_id);
    void setValue_BoundaryVolumeNode(const uiint& mesh_id, const uiint& bnd_id, const uiint& bnode_id, const uiint& dof, const double& val);
    void GeneBoundaryVolume(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType, const uiint& elemType,
                            const uiint& mesh_id, const uiint& elem_id, vuint& vBNodeID,
                            const uiint& b_vol_id, const uiint& dof, const double& val);
    void initEdgeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void resizeEdgeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void GeneBoundaryEdgeNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
                              const uiint& mesh_id, const uiint& node_id,
                              const uiint& b_node_id);
    void setValue_BoundaryEdgeNode(const uiint& mesh_id, const uiint& bnd_id, const uiint& bnode_id, const uiint& dof, const double& val);
    void GeneBoundaryEdge(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType, const uiint& elemType,
                          const uiint& mesh_id, const uiint& elem_id, const uiint& edge_id, vuint& vBNodeID,
                          const uiint& b_edge_id, const uiint& dof, const double& val);

    void refineBoundary();

    void initBucketNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& maxID, const uiint& minID);
    void setIDBucketNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const uiint& index);
    void setupBucketNode(const uiint& mgLevel, const uiint& mesh_id);
    void initBucketElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& maxID, const uiint& minID);
    void setIDBucketElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const uiint& index);
    void setupBucketElement(const uiint& mgLevel, const uiint& mesh_id);
    void setupBucketMesh(const uiint& mgLevel, const uiint& maxID, const uiint& minID);

    void reserveMaterial(const uiint& res_size);
    void GeneMaterial(const uiint& mesh_id, const uiint& material_id, string& name, vuint& vType, vdouble& vValue);
    void reserveCommMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& res_size);
    void GeneCommMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& comID, const uiint& myRank, const uiint& nTransmitRank);
    void reserveCommNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& commesh_id, const uiint& res_size);
    void GeneCommNode(const uiint& mgLevel, const uiint& commNodeID,
                      const uiint& mesh_id, const uiint& commesh_id, const uiint& nodeID, const uiint& rank);
    void reserveCommElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& commesh_id, const uiint& res_size);
    void GeneCommElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& commesh_id,
                         const uiint& nType, const uiint& elemID, vuint& vCommNodeID);
protected:
    void GeneProgCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem);
    void dividCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem);
public:
    //--
    // ContactMesh
    //--
    void GeneContactMesh(const uiint& contactID, const uiint& myRank, const uiint& nProp);
    void GeneContactNode(const uiint& mgLevel, const uiint& contactID, const uiint& conNodeID, const vdouble& vCoord,
                         const string& s_param_type, const uiint& numOfVector, const uiint& numOfScalar,
                         bool bmesh, const uiint& meshID, const uiint& nodeID,
                         const uiint& rank, const uiint& maslave);
    void GeneMasterFace(const uiint& contactID, const uiint& shapeType, const uiint& masterFaceID,
                        bool bmesh, const uiint& meshID, const uiint& elemID, const uiint& elemFaceID,
                        const vuint& vConNodeID, const uiint& face_rank);
    void GeneSlaveFace(const uiint& contactID, const uiint& shapeType, const uiint& slaveFaceID,
                       bool bmesh, const uiint& meshID, const uiint& elemID, const uiint& elemFaceID,
                       const vuint& vConNodeID, const uiint& face_rank);
    void refineContactMesh();
protected:
    void setProgHexaMPCMaster(vector<CElement*>& vProgElem, const uiint& iface, const uiint& i, const uiint& j, const uiint& k, const uiint& l);
    void setProgHexaMPCSlave(vector<CElement*>& vProgElem, const uiint& iface, const uiint& i, const uiint& j, const uiint& k, const uiint& l);
    // CommMesh2属性設定(面 辺 点)
    void setProgHexaCommMesh2Face(vector<CElement*>& vProgElem, const uiint& iface, const uiint& i, const uiint& j, const uiint& k, const uiint& l);
    void setProgHexaCommMesh2Edge(vector<CElement*>& vProgElem, const uiint& iedge, const uiint& i, const uiint& j);
    void setProgHexaCommMesh2Vert(vector<CElement*>& vProgElem, const uiint& ivert);
public:
    void refineCommMesh2();
    void GeneCommMesh2(const uiint& mgLevel, const uiint& mesh_id, const uiint& comID, const uiint& numOfFace, const uiint& numOfCommNode,
                       const uiint& myRank, const uiint& nTransmitRank);
    void GeneCommFace(const uiint& mgLevel, const uiint& commeshID, const uiint& face_id,
                      const uiint& mesh_id,const uiint elem_id, const uiint& elem_ent_num, const uiint& elem_type, const vuint& vCommNodeID);
    void GeneCommNodeCM2(const uiint& mgLevel, const uiint& mesh_id, const uiint& node_id, const uiint& commeshID, const uiint& comm_node_id, const vdouble& vCoord);

    void GeneElemGrpOBJ(const uiint& mgLevel, const uiint& mesh_id, const vuint& vGrpID, vstring& vGrpName);
    void setElemID_with_ElemGrp(const uiint& mgLevel, const uiint& mesh_id, const uiint& nGrpID, const vuint& vElemID);

    void setupBNodeMarking();//----------------- 境界Nodeマーキング for CMesh : 配列の生成&マーキング両方を行う
    void setupLargeRankCommNode_Marking();//---- Rank大の通信ノードのマーキング for CMesh : 配列の生成&マーキング両方を行う
};
#endif
}
