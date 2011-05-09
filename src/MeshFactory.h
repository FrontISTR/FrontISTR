//
//  MeshFactory.h
//
//			2009.07.23
//			2008.11.05
//			k.Takeda
#ifndef MESH_FACTORY_HH_E74CAEBF_35A3_407b_A7B8_49D7A2DDDF8E
#define MESH_FACTORY_HH_E74CAEBF_35A3_407b_A7B8_49D7A2DDDF8E

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
//#include "IndexBucketMesh.h"

#include "Material.h"

#include "CommunicationMesh.h"// CommElementによる,通信領域Mesh
#include "ProgElementTree.h"// CommElementのprolongation時に利用

#include "boost/lexical_cast.hpp"//数=>文字変換

// MPC
#include "SkinFace.h"
#include "MasterFace.h"
#include "ContactMesh.h"
#include "ContactNode.h"

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
    //CMeshFactory(void);
    virtual ~CMeshFactory(void);

protected:
    CGMGModel    *mpGMGModel;

    CAssyModel   *mpTAssyModel;// 各 AssembledMesh Model
    CMesh        *mpTMesh;     // 各 Mesh
    CCommMesh    *mpTCommMesh; // 各 CommMesh(通信領域) in Mesh

    uint    mMGLevel;//マルチグリッド階層数 => 1.ならマルチグリッドなし

    Utility::CLogger *mpLogger;

public:
    // root Model
    virtual void setGMGModel(CGMGModel* pGMGModel){ mpGMGModel = pGMGModel;}
    // assembled Model
    void GeneAssyModel(const uint& num_of_mgLevel);// AssyModelをmgLevel数ぶんだけ生成
    // each Mesh in AssyModel at mgLevel
    void reserveMesh(const uint& mgLevel,const uint& num_of_mesh);
    void GeneMesh(const uint& mgLevel, const uint& mesh_id, const uint& index);// MeshをAssyModelに生成 at mgLevel


    // Refine for MultiGrid
    // --------------------
    void setMGLevel(const uint& num_of_mglevel){ mMGLevel= num_of_mglevel;}
    uint& getMGLevel(){ return mMGLevel;}
    void refineMesh();// <<<<<<<<<<<<< -- MultiGrid prolongation.
protected:
    void GeneProgElem(const uint& ilevel,CElement* pElem, CElement* pProgElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh);
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
    void reserveAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node);
    void GeneAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node);


    // Boundary
    // --------
    void reserveBoundaryNode(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void reserveBoundaryFace(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void reserveBoundaryVolume(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd);
    void GeneBoundaryNode(const uint& mgLevel, const uint& mesh_id, const uint& id,
                          const uint& dof, const uint& bndType,const vdouble& vVal);
    void GeneBoundaryFace(const uint& mgLevel, const uint& mesh_id, const uint& elem_id, const uint& face_id,
                          const uint& dof, const uint& bndType,const vdouble& vVal);
    void GeneBoundaryVolume(const uint& mgLevel, const uint& mesh_id, const uint& id,
                          const uint& dof, const uint& bndType,const vdouble& vVal);



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
    void GeneContactMesh(const uint& contactID);//接合メッシュの生成:ContactMesh数ぶんコール
    void GeneContactNode(const uint& mgLevel, const uint& contactID, const uint& conNodeID, const vdouble& vCoord,
            const string& s_param_type, const uint& numOfVector, const uint& numOfScalar,
            bool bmesh, const uint& meshID, const uint& nodeID,
            const uint& rank, const uint& maslave);
    void GeneMasterFace(const uint& contactID, const uint& shapeType, const uint& masterFaceID,
            bool bmesh, const uint& meshID, const uint& elemID, const uint& elemFaceID,
            const vuint& vConNodeID);
    void GeneSlaveFace(const uint& contactID, const uint& shapeType, const uint& slaveFaceID,
            bool bmesh, const uint& meshID, const uint& elemID, const uint& elemFaceID,
            const vuint& vConNodeID);

    void refineContactMesh();

protected:
    void setProgHexaMPCMaster(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l);
    void setProgHexaMPCSlave(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l);
    
};
}
#endif
