/* 
 * File:   BoundaryVolumeMesh.h
 * Author: ktakeda
 *
 * 要素境界(BoundaryVolume)管理
 *
 * Created on 2010/04/12, 17:00
 */
#include "TypeDef.h"
#include "BoundaryVolume.h"//個別の体積境界(値)
#include "BoundaryMesh.h"

#include "ElementType.h"

// Equivalent Nodal Force(等価節点力)
#include "ShapeHexa.h"
#include "ShapeTetra.h"
#include "ShapePrism.h"

namespace pmw{
#ifndef _BOUNDARYVOLUMEMESH_H
#define	_BOUNDARYVOLUMEMESH_H
class CBoundaryVolumeMesh:public CBoundaryMesh{
public:
    CBoundaryVolumeMesh();
    virtual ~CBoundaryVolumeMesh();

protected:
    vector<CBoundaryVolume*> mvBVolume;
    map<uint, uint, less<uint> > mmBVolumeID2Index;

    //(辺-面-体)にBNodeを生成して、Refineで使用
    vector<CBoundaryNode*> mvBEdgeBNode;
    vector<CBoundaryNode*> mvBFaceBNode;
    vector<CBoundaryNode*> mvBVolBNode;

public:
    //境界体積 :=> BVolume
    void resizeVolume(const uint& res_size);
    uint getNumOfVolume(){ return mvBVolume.size();}
    
    void setBVolume(const uint& index, CBoundaryVolume *pBVolume);
    void addBVolume(CBoundaryVolume *pBVolume);

    CBoundaryVolume* getBVolumeIX(const uint& index){ return mvBVolume[index];}
    CBoundaryVolume* getBVolumeID(const uint& id){
        uint index= mmBVolumeID2Index[id];
        return mvBVolume[index];
    }
    uint& getBVolumeIndex(const uint& id){ mmBVolumeID2Index[id];}


protected:
    vvuint mvAggregateVol;//頂点のVol-IDの集合

public:
    void resizeAggVol();
    void setupAggVol();
    void setAggVol(const uint& ibnode, const uint& nVolumeID){ mvAggregateVol[ibnode].push_back(nVolumeID);}

    void GeneEdgeBNode();
    void GeneFaceBNode();
    void GeneVolBNode();

//    CBoundaryNode* getEdgeBNode(const uint& iedge){  return mvBEdgeBNode[iedge];}// 辺BNode提供
//    CBoundaryNode* getFaceBNode(const uint& iface){  return mvBFaceBNode[iface];}// 面BNode提供
//    CBoundaryNode* getVolumeBNode(const uint& ivol){ return mvBVolBNode[ivol];}  // 体BNode提供

    // 境界要素の再分割 => progBMeshにセット
    // ----
    void refine(CBoundaryVolumeMesh *pProgBVolMesh);
    
    void deleteProgData();// Refine 後処理 : 辺-面 BNode vectorの解放

protected:
    //Neumann条件の節点分配(形状関数による等価節点力:EquivalentNodeForce)
    //----
    virtual void distNeumannValue();

    //Dirichlet条件の節点分配
    // Level==0:BNode周囲の要素集合平均値
    // Level>=1:BNode間の平均
    //----
    virtual void distDirichletValue_at_CGrid();//Level==0(coarse grid)
    virtual void distDirichletValue_at_FGrid();//Level>=1(fine grid)
};
#endif	/* _BOUNDARYVOLUMEMESH_H */
}





