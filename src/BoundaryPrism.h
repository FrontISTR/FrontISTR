/* 
 * File:   BoundaryPrism.h
 * Author: ktakeda
 *
 * Created on 2010/07/02, 17:10
 */
#include "BoundaryVolume.h"
#include "BoundaryHexa.h"//refineで使用

namespace pmw{
#ifndef _BOUNDARYPRISM_H
#define	_BOUNDARYPRISM_H
class CBoundaryPrism:public CBoundaryVolume{
public:
    CBoundaryPrism();// 形状( Tetra | Hexa | Prism )初期化
    virtual ~CBoundaryPrism();
private:
    static uint mnElemType;
    static uint mNumOfFace;
    static uint mNumOfEdge;
    static uint mNumOfNode;

protected:
    virtual uint* getLocalNode_Edge(const uint& iedge);
    virtual uint* getLocalNode_Face(const uint& iface);

public:
    virtual uint getElemType();
    virtual uint getNumOfEdge();
    virtual uint getNumOfFace();

    
    virtual PairBNode getPairBNode(const uint& iedge);
    virtual uint& getEdgeID(PairBNode& pairBNode);

    virtual vector<CBoundaryNode*> getFaceCnvNodes(const uint& iface);
    virtual uint& getFaceID(vector<CBoundaryNode*>& vBNode);
    
    
    virtual void refine(uint& countID, const vuint& vDOF);// Refine 再分割
    virtual double& calcVolume();// BoundaryVolumeの体積

    virtual void distDirichletVal(const uint& dof, const uint& mgLevel);//上位グリッドBNodeへのディレクレ値の分配

    virtual void deleteProgData();// Refine 後処理 : 辺-面 BNode vectorの解放
};
#endif	/* _BOUNDARYPRISM_H */
}









