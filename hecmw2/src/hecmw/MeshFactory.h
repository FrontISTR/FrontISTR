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
    

    uiint    mMGLevel;//マルチグリッド階層数 => 1.ならマルチグリッドなし
    uiint mnSolutionType;//FEM, FVM

    Utility::CLogger *mpLogger;

public:
    // root Model
    virtual void setGMGModel(CGMGModel* pGMGModel){ mpGMGModel = pGMGModel;}
    // assembled Model
    void GeneAssyModel(const uiint& nNumOfLevel);// AssyModelをmgLevel数ぶんだけ生成
    // each Mesh in AssyModel at mgLevel
    void reserveMesh(const uiint& mgLevel,const uiint& num_of_mesh);
    void GeneMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& index, const uiint& nProp);// MeshをAssyModelに生成 at mgLevel


    // Solution Type :: FEM .or. FVM
    void setSolutionType(const uiint& nSolutionType){ mnSolutionType = nSolutionType;}

    
    // Refine for MultiGrid
    // --------------------
    void setMGLevel(const uiint& num_of_mglevel){ mMGLevel= num_of_mglevel;}
    uiint& getMGLevel(){ return mMGLevel;}
    void refineMesh();// <<<<<<<<<<<<< -- MultiGrid prolongation.
protected:
    void MGMeshConstruct();// <<<<<<<<<<<<< -- MultiGrid prolongation.
    void SGMeshConstruct();// <<<<<<<<<<<<< -- SingleGrid
    void GeneProgElem(const uiint& ilevel,CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividHexa(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividTetra(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividPrism(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividPyramid(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividQuad(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividTriangle(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);
    void dividBeam(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh);

public:
    // Node
    // ----
    void reserveNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_node);
    void GeneNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const vdouble& coord,
                  const uiint& nodeType, const uiint& nNumOfSDOF, const uiint& nNumOfVDOF);
    void setupNode(const uiint& mgLevel, const uiint& mesh_id);//numOfNodeのセット


    // Element
    // -------
    void reserveElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_element);
    void GeneElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const uiint& type, const vint& node_id);
    void setupElement(const uiint& mgLevel, const uiint& mesh_id);//numOfElementのセット
    

    // Aggregate-Element, Aggregate-Node
    //
    void resizeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_node);
    void GeneAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_node);


    // BoundaryMesh
    // --------
    void reserveBoundaryNodeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd);
    void GeneBoundaryNodeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name);
    void reserveBoundaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd);
    void GeneBoundaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name,
                              const uiint& numOfDOF, const vuint& vDOF);
    void reserveBoundaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd);
    void GeneBoundaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name,
                                const uiint& numOfDOF, const vuint& vDOF);
    void reserveBoundaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd);
    void GeneBoundaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name,
                              const uiint& numOfDOF, const vuint& vDOF);
    uiint getNumOfBounaryNodeMesh(const uiint& mgLevel, const uiint& mesh_id);
    uiint getNumOfBounaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id);
    uiint getNumOfBounaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id);
    uiint getNumOfBounaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id);
    CBoundaryFaceMesh* getBoundaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    CBoundaryVolumeMesh* getBoundaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    CBoundaryEdgeMesh* getBoundaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);

    // 境界節点
    void GeneBoundaryNode(const uiint& mgLevel, const uiint& bnd_id,  const uiint& bndType,
                          const uiint& mesh_id, const uiint& node_id,
                          const uiint& b_node_id, const uiint& dof, const double& val);
    // 境界面メッシュ-節点
    void initFaceAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void resizeFaceAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void GeneBoundaryFaceNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
                          const uiint& mesh_id, const uiint& node_id,
                          const uiint& b_node_id);
    void setValue_BoundaryFaceNode(const uiint& mesh_id, const uiint& bnd_id, const uiint& bnode_id, const uiint& dof, const double& val);//新Dirichlet用途
    // 境界面メッシュ-面
    void GeneBoundaryFace(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType, const uiint& elemType,
                          const uiint& mesh_id, const uiint& elem_id, const uiint& face_id, vuint& vBNodeID,
                          const uiint& b_face_id, const uiint& dof, const double& val);
    // 境界体積メッシュ-節点
    void initVolumeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void resizeVolumeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void GeneBoundaryVolumeNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
                          const uiint& mesh_id, const uiint& node_id,
                          const uiint& b_node_id);
    void setValue_BoundaryVolumeNode(const uiint& mesh_id, const uiint& bnd_id, const uiint& bnode_id, const uiint& dof, const double& val);//新Dirichlet用途
    // 境界体積メッシュ-体積
    void GeneBoundaryVolume(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType, const uiint& elemType,
                            const uiint& mesh_id, const uiint& elem_id, vuint& vBNodeID,
                            const uiint& b_vol_id, const uiint& dof, const double& val);
    // 境界辺メッシュ-節点
    void initEdgeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void resizeEdgeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id);
    void GeneBoundaryEdgeNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
                          const uiint& mesh_id, const uiint& node_id,
                          const uiint& b_node_id);
    void setValue_BoundaryEdgeNode(const uiint& mesh_id, const uiint& bnd_id, const uiint& bnode_id, const uiint& dof, const double& val);//新Dirichlet用途
    // 境界辺メッシュ-辺
    void GeneBoundaryEdge(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType, const uiint& elemType,
                          const uiint& mesh_id, const uiint& elem_id, const uiint& edge_id, vuint& vBNodeID,
                          const uiint& b_edge_id, const uiint& dof, const double& val);
    // ----
    // 境界条件階層化
    // ----
    void refineBoundary();



    // Bucket(in Mesh) setup
    // ---------------------
    // Node
    void initBucketNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& maxID, const uiint& minID);
    void setIDBucketNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const uiint& index);
    void setupBucketNode(const uiint& mgLevel, const uiint& mesh_id);// Nodeオブジェクト設定後に一度にBucketをセットアップ
    // Element
    void initBucketElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& maxID, const uiint& minID);
    void setIDBucketElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const uiint& index);
    void setupBucketElement(const uiint& mgLevel, const uiint& mesh_id);// Elementオブジェクト設定後に一度にBucketをセットアップ

    // Bucket(in AssyModel) setup
    // --------------------------
    void setupBucketMesh(const uiint& mgLevel, const uiint& maxID, const uiint& minID);


    // Material
    // --
    void reserveMaterial(const uiint& res_size);
    void GeneMaterial(const uiint& mesh_id, const uiint& material_id, string& name, vuint& vType, vdouble& vValue);


    // 通信領域(CommunicationMesh)
    // --
    void reserveCommMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& res_size);
    void GeneCommMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& comID, const uiint& myRank, const uiint& nTransmitRank);
    // CommNode
    void reserveCommNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& commesh_id, const uiint& res_size);
    void GeneCommNode(const uiint& mgLevel, const uiint& commNodeID,
                      const uiint& mesh_id, const uiint& commesh_id, const uiint& nodeID, const uiint& rank);
    // CommElement
    void reserveCommElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& commesh_id, const uiint& res_size);
    void GeneCommElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& commesh_id,
                                   const uiint& nType, const uiint& elemID, vuint& vCommNodeID);
    
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
    void GeneContactMesh(const uiint& contactID, const uiint& myRank, const uiint& tranRank, const uiint& nProp);//接合メッシュの生成:ContactMesh数ぶんコール
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
    // MPC属性をprogElemにセット :Hexa専用
    void setProgHexaMPCMaster(vector<CElement*>& vProgElem, const uiint& iface, const uiint& i, const uiint& j, const uiint& k, const uiint& l);
    void setProgHexaMPCSlave(vector<CElement*>& vProgElem, const uiint& iface, const uiint& i, const uiint& j, const uiint& k, const uiint& l);

    // 通信界面(要素分割:節点共有型)
    //
    // CommMesh2(節点-通信界面)の属性をprogElemにセット :Hexa専用
    void setProgHexaCommMesh2Entity(vector<CElement*>& vProgElem, const uiint& iface, const uiint& i, const uiint& j, const uiint& k, const uiint& l);
public:
    void refineCommMesh2();
    void GeneCommMesh2(const uiint& mgLevel, const uiint& mesh_id, const uiint& comID, const uiint& numOfFace, const uiint& numOfCommNode,
            const uiint& myRank, const uiint& nTransmitRank);
    void GeneCommFace(const uiint& mgLevel, const uiint& commeshID, const uiint& face_id,
            const uiint& mesh_id,const uiint elem_id, const uiint& elem_ent_num, const uiint& elem_type, const vuint& vCommNodeID);
    void GeneCommNodeCM2(const uiint& mgLevel, const uiint& mesh_id, const uiint& node_id, const uiint& commeshID, const uiint& comm_node_id, const vdouble& vCoord);//CommMesh2用途


    // --
    // グループ
    // --
    // ・GroupObjectの生成
    // ・GroupID, GroupNameのセット
    void GeneElemGrpOBJ(const uiint& mgLevel, const uiint& mesh_id, const vuint& vGrpID, vstring& vGrpName);//ElementGroupの生成
    // ・指定GrpIDへパラメーターをセット
    void setElemID_with_ElemGrp(const uiint& mgLevel, const uiint& mesh_id, const uiint& nGrpID, const vuint& vElemID);


    ////    // --
    ////    // Resデータ(Nodeへのセット) :=> 使用中止
    ////    // --
    ////    void setNodeValue(const uiint& mgLevel, const uiint& nMeshID, const uiint& nNodeID,
    ////                        const uiint& nNumOfSDOF, const uiint& nNumOfVDOF, vdouble& vScaValue, vdouble& vVecValue);
    
};
#endif  //_MESH_FACTORY_HH_E74C
}

