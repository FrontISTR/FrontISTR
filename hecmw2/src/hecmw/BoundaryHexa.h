/* 
 * File:   BoundaryHexa.h
 * Author: ktakeda
 *
 * Created on 2010/07/02, 17:10
 */
#include "BoundaryVolume.h"

namespace pmw{
#ifndef _BOUNDARYHEXA_H
#define	_BOUNDARYHEXA_H
class CBoundaryHexa:public CBoundaryVolume{
public:
    CBoundaryHexa();// 形状( Tetra | Hexa | Prism )初期化
    virtual ~CBoundaryHexa();
//private:
//    static uint mnElemType;
//    static uint mNumOfFace;
//    static uint mNumOfEdge;
//    static uint mNumOfNode;

protected:
    virtual uiint* getLocalNode_Edge(const uiint& iedge);
    virtual uiint* getLocalNode_Face(const uiint& iface);

public:
    virtual uiint getElemType();
    virtual uiint getNumOfEdge();
    virtual uiint getNumOfFace();
    virtual uiint getNumOfNode();
    virtual uiint getNumOfVert();

    virtual void setOrder(const uiint& order);

    virtual PairBNode getPairBNode(const uiint& iedge);
    virtual uiint& getEdgeID(PairBNode& pairBNode);

    virtual vector<CBoundaryNode*> getFaceCnvNodes(const uiint& iface);
    virtual uiint& getFaceID(vector<CBoundaryNode*>& vBNode);



    virtual void refine(uiint& countID, const vuint& vDOF);// Refine 再分割

    virtual double& calcVolume();// BoundaryVolumeの体積

    virtual void distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel);//上位グリッドBNodeへのディレクレ値の分配

    virtual void replaceEdgeBNode(const uiint& iedge);//2次要素の場合に辺BNodeをmvBNodeへ移設.
    
    virtual void deleteProgData();// Refine 後処理 : 辺-面 BNode vectorの解放
};
#endif	/* _BOUNDARYHEXA_H */
}




