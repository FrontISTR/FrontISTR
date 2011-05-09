//
//  MeshFactory.h
//
//			2009.07.23
//			2008.11.05
//			k.Takeda
#include "CommonStd.h"
#include "TypeDef.h"

#include "Logger.h"
#include "Element4Factory.h"// Elemnt => Hexa,Tetra,Prism,Pyramid,Quad,Triangle,Beam
                            // CommElement => CommHexa, CommTetra, CommPrism
#include "ElementProperty.h"// 頂点数,面数,辺数

#include "ScalarNode.h"
#include "VectorNode.h"
#include "ScalarVectorNode.h"
#include "NodeType.h"


// 境界Mesh
#include "BoundaryNodeMesh.h"
#include "BoundaryFaceMesh.h"
#include "BoundaryVolumeMesh.h"
#include "BoundaryEdgeMesh.h"
// 境界オブジェクト
#include "BoundaryNode.h"
#include "BoundarySBNode.h"//BoundaryNodeMesh用途
#include "BoundaryFace.h"
#include "BoundaryHexa.h"
#include "BoundaryTetra.h"
#include "BoundaryPrism.h"
#include "BoundaryEdge.h"
#include "BoundaryType.h"
// 境界BNodeMesh
#include "BNodeMeshGrp.h"


#include "AggregateElement.h"
#include "AggregateNode.h"
#include "Mesh.h"
#include "AssyModel.h"
#include "GMGModel.h"

#include "IndexBucket.h"
//#include "IndexBucketMesh.h"

#include "Material.h"

// CommMesh(要素共有型)
#include "CommunicationMesh.h"// CommElementによる,通信領域Mesh
#include "ProgElementTree.h"// CommElementのprolongation時に利用

#include "boost/lexical_cast.hpp"//数=>文字変換

// MPC
#include "SkinFace.h"
#include "MasterFace.h"
#include "ContactMesh.h"
#include "ContactNode.h"


// CommMesh2(節点共有型)
#include "CommMesh2.h"
#include "CommFace.h"
#include "CommNode.h"

// SolutionType
#include "SolutionType.h"


// グループ
#include "ElementGroup.h"

//MPI
#include "HEC_MPI.h"

namespace pmw{
#ifndef _MESH_FACTORY_HH_E74C
#define _MESH_FACTORY_HH_E74C
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
    //CMeshFactory(void);
    virtual ~CMeshFactory(void);

protected:
    CGMGModel    *mpGMGModel;

    CAssyModel   *mpTAssyModel;// 各 AssembledMesh Model
    CMesh        *mpTMesh;     // 各 Mesh
    CCommMesh    *mpTCommMesh; // 各 CommMesh (要素共有型-通信) in Mesh
    CCommMesh2   *mpTCommMesh2;// 各 CommMesh2(節点共有型-通信) in Mesh
    

    uint    mMGLevel;//マルチグリッド階層数 => 1.ならマルチグリッドなし
    uint mnSolutionType;//FEM, FVM

    Utility::CLogger *mpLogger;

public:
    // root Model
    virtual void setGMGModel(CGMGModel* pGMGModel){ mpGMGModel = pGMGModel;}
    // assembled Model
    void GeneAssyModel(const uint& num_of_mgLevel);// AssyModelをmgLevel数ぶんだけ生成
    // each Mesh in AssyModel at mgLevel
    void reserveMesh(const uint& mgLevel,const uint& num_of_mesh);
    void GeneMesh(const uint& mgLevel, const uint& mesh_id, const uint& index, const uint& nProp);// MeshをAssyModelに生成 at mgLevel


    // Solution Type :: FEM .or. FVM
    void setSolutionType(const uint& nSolutionType){ mnSolutionType = nSolutionType;}

    
    // Refine for MultiGrid
    // --------------------
    void setMGLevel(const uint& num_of_mglevel){ mMGLevel= num_of_mglevel;}
    uint& getMGLevel(){ return mMGLevel;}
    void refineMesh();// <<<<<<<<<<<<< -- MultiGrid prolongation.
protected:
    void MGMeshConstruct();// <<<<<<<<<<<<< -- MultiGrid prolongation.
    void SGMeshConstruct();// <<<<<<<<<<<<< -- SingleGrid
    void GeneProgElem(const uint& ilevel,CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividHexa(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividTetra(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividPrism(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividPyramid(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividQuad(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividTriangle(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
    void dividBeam(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);

public:
    // Node
    // ----
    void reserveNode(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node);
    void GeneNode(const uint& mgLevel, const uint& mesh_id, const uint& id, const vdouble& coord,
                  const uint& nodeType, const uint& numOfScaParam, const uint& numOfVecParam);
    void setupNode(const uint& mgLevel, const uint& mesh_id);//numOfNodeのセット


    // Element
    // -------
    void reserveElement(const uint& mgLevel, const uint& mesh_id, const uint& num_of_element);
    void GeneElement(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& type, const vint& node_id);
    void setupElement(const uint& mgLevel, const uint& mesh_id);//numOfElementのセット    
    

    // Aggregate-Element, Aggregate-Node
    //
    void resizeAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node);
    void GeneAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node);


    // BoundaryMesh
    // --------
    void reserveBoundaryNodeMesh(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void GeneBoundaryNodeMesh(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id, const uint& bnd_type, const string& bnd_name);
    void reserveBoundaryFaceMesh(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void GeneBoundaryFaceMesh(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id, const uint& bnd_type, const string& bnd_name,
                              const uint& numOfDOF, const vuint& vDOF);
    void reserveBoundaryVolumeMesh(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void GeneBoundaryVolumeMesh(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id, const uint& bnd_type, const string& bnd_name,
                                const uint& numOfDOF, const vuint& vDOF);
    void reserveBoundaryEdgeMesh(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void GeneBoundaryEdgeMesh(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id, const uint& bnd_type, const string& bnd_name,
                              const uint& numOfDOF, const vuint& vDOF);
    uint getNumOfBounaryNodeMesh(const uint& mgLevel, const uint& mesh_id);
    uint getNumOfBounaryFaceMesh(const uint& mgLevel, const uint& mesh_id);
    uint getNumOfBounaryVolumeMesh(const uint& mgLevel, const uint& mesh_id);
    uint getNumOfBounaryEdgeMesh(const uint& mgLevel, const uint& mesh_id);

    // 境界節点
    void GeneBoundaryNode(const uint& mgLevel, const uint& bnd_id,  const uint& bndType,
                          const uint& mesh_id, const uint& node_id,
                          const uint& b_node_id, const uint& dof, const double& val);
    // 境界面メッシュ-節点
    void initFaceAggregate(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id);
    void GeneBoundaryFaceNode(const uint& mgLevel, const uint& bnd_id, const uint& bndType,
                          const uint& mesh_id, const uint& node_id,
                          const uint& b_node_id);
    // 境界面メッシュ-面
    void GeneBoundaryFace(const uint& mgLevel, const uint& bnd_id, const uint& bndType, const uint& elemType,
                          const uint& mesh_id, const uint& elem_id, const uint& face_id, vuint& vBNodeID,
                          const uint& b_face_id, const uint& dof, const double& val);
    // 境界体積メッシュ-節点
    void initVolumeAggregate(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id);
    void GeneBoundaryVolumeNode(const uint& mgLevel, const uint& bnd_id, const uint& bndType,
                          const uint& mesh_id, const uint& node_id,
                          const uint& b_node_id);
    // 境界体積メッシュ-体積
    void GeneBoundaryVolume(const uint& mgLevel, const uint& bnd_id, const uint& bndType, const uint& elemType,
                            const uint& mesh_id, const uint& elem_id, vuint& vBNodeID,
                            const uint& b_vol_id, const uint& dof, const double& val);
    // 境界辺メッシュ-節点
    void initEdgeAggregate(const uint& mgLevel, const uint& mesh_id, const uint& bnd_id);
    void GeneBoundaryEdgeNode(const uint& mgLevel, const uint& bnd_id, const uint& bndType,
                          const uint& mesh_id, const uint& node_id,
                          const uint& b_node_id);
    // 境界辺メッシュ-辺
    void GeneBoundaryEdge(const uint& mgLevel, const uint& bnd_id, const uint& bndType, const uint& elemType,
                          const uint& mesh_id, const uint& elem_id, const uint& edge_id, vuint& vBNodeID,
                          const uint& b_edge_id, const uint& dof, const double& val);
    // ----
    // 境界条件階層化
    // ----
    void refineBoundary();



    // Bucket(in Mesh) setup
    // ---------------------
    // Node
    void initBucketNode(const uint& mgLevel, const uint& mesh_id, const uint& maxID, const uint& minID);
    void setIDBucketNode(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& index);
    void setupBucketNode(const uint& mgLevel, const uint& mesh_id);// Nodeオブジェクト設定後に一度にBucketをセットアップ
    // Element
    void initBucketElement(const uint& mgLevel, const uint& mesh_id, const uint& maxID, const uint& minID);
    void setIDBucketElement(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& index);
    void setupBucketElement(const uint& mgLevel, const uint& mesh_id);// Elementオブジェクト設定後に一度にBucketをセットアップ

    // Bucket(in AssyModel) setup
    // --------------------------
    void setupBucketMesh(const uint& mgLevel, const uint& maxID, const uint& minID);


    // Material
    // --
    void reserveMaterial(const uint& res_size);
    void GeneMaterial(const uint& mesh_id, const uint& material_id, string& name, vuint& vType, vdouble& vValue);


    // 通信領域(CommunicationMesh)
    // --
    void reserveCommMesh(const uint& mgLevel, const uint& mesh_id, const uint& res_size);
    void GeneCommMesh(const uint& mgLevel, const uint& mesh_id, const uint& comID, const uint& myRank, const uint& nTransmitRank);
    // CommNode
    void reserveCommNode(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, const uint& res_size);
    void GeneCommNode(const uint& mgLevel, const uint& commNodeID, 
                      const uint& mesh_id, const uint& commesh_id, const uint& nodeID, const uint& rank);
    // CommElement
    void reserveCommElement(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, const uint& res_size);
    void GeneCommElement(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, 
                                   const uint& nType, const uint& elemID, vuint& vCommNodeID);
    
protected:
    // CommMeshのprolongation: Meshのprolongation処理に組み込み
    // --
    // refineMeshからコールされる,CommElemの再分割処理部分// <<<<<<<<< MultiGrid prolongation for CommMesh
    // --
    void GeneProgCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem);//1個のCommElemに生成するprogCommElemの生成
    void dividCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem);//CommElemの分割
    
public:
    // MPCの為の表面メッシュ(Skin)生成
    // --
    void GeneContactMesh(const uint& contactID, const uint& myRank, const uint& tranRank, const uint& nProp);//接合メッシュの生成:ContactMesh数ぶんコール
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
    // MPC属性をprogElemにセット :Hexa専用
    void setProgHexaMPCMaster(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l);
    void setProgHexaMPCSlave(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l);

    // 通信界面(要素分割:節点共有型)
    //
    // CommMesh2(節点-通信界面)の属性をprogElemにセット :Hexa専用
    void setProgHexaCommMesh2Entity(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l);
public:
    void refineCommMesh2();
    void GeneCommMesh2(const uint& mgLevel, const uint& mesh_id, const uint& comID, const uint& numOfFace, const uint& numOfCommNode,
            const uint& myRank, const uint& nTransmitRank);
    void GeneCommFace(const uint& mgLevel, const uint& commeshID, const uint& face_id,
            const uint& mesh_id,const uint elem_id, const uint& elem_ent_num, const uint& elem_type, const vuint& vCommNodeID);
    void GeneCommNodeCM2(const uint& mgLevel, const uint& mesh_id, const uint& node_id, const uint& commeshID, const uint& comm_node_id, const vdouble& vCoord);//CommMesh2用途


    // --
    // グループ
    // --
    // ・GroupObjectの生成
    // ・GroupID, GroupNameのセット
    void GeneElemGrpOBJ(const uint& mgLevel, const uint& mesh_id, const vuint& vGrpID, vstring& vGrpName);//ElementGroupの生成
    // ・指定GrpIDへパラメーターをセット
    void setElemID_with_ElemGrp(const uint& mgLevel, const uint& mesh_id, const uint& nGrpID, const vuint& vElemID);

};
#endif  //_MESH_FACTORY_HH_E74C
}

