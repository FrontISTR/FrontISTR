/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ProgElementTree.h
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
#include "TypeDef.h"
#include "ElementType.h"
#include "Logger.h"
namespace pmw
{
#ifndef _PROGELEMENTTREE_H
#define	_PROGELEMENTTREE_H
class CProgElementTree
{
private:
    CProgElementTree();
public:
    static CProgElementTree* Instance() {
        static CProgElementTree progElemTree;
        return &progElemTree;
    }
    virtual ~CProgElementTree();
protected:
    uiint mInvalidNum;
    uiint mProgVert;
    uiint mHexaVertChildVert[8];
    uiint mTetraVertChildVert[4];
    uiint mPrismVertChildVert[6];
    uiint mPyramidVertChildVert[5];
    uiint mQuadVertChildVert[4];
    uiint mTriangleVertChildVert[3];
    uiint mBeamVertChildVert[2];
    uiint mHexaEdgeChildVert[8][12];
    uiint mTetraEdgeChildVert[4][6];
    uiint mPrismEdgeChildVert[6][9];
    uiint mPyramidEdgeChildVert[8][8];
    uiint mQuadEdgeChildVert[4][4];
    uiint mTriangleEdgeChildVert[3][3];
    uiint mBeamEdgeChildVert[2][1];
    uiint mHexaFaceChildVert[8][6];
    uiint mTetraFaceChildVert[4][4];
    uiint mPrismFaceChildVert[6][5];
    uiint mPyramidFaceChildVert[8][5];
    uiint mQuadFaceChildVert[4][1];
    uiint mTriangleFaceChildVert[3][1];
    uiint mHexaVolChildVert[8];
    uiint mTetraVolChildVert[4];
    uiint mPrismVolChildVert[6];
    uiint mPyramidVolChildVert[8];
    uiint mQuadVolChildVert[4];
    uiint mTriangleVolChildVert[3];
    uiint mBeamVolChildVert[2];
public:
    uiint& getInvalidNum() {
        return mInvalidNum;
    }
    uiint& getHexaVertProgVert(const uiint& child_address) {
        return mHexaVertChildVert[child_address];
    }
    uiint& getTetraVertProgVert(const uiint& child_address) {
        return mTetraVertChildVert[child_address];
    }
    uiint& getPrismVertProgVert(const uiint& child_address) {
        return mPrismVertChildVert[child_address];
    }
    uiint& getPyramidVertProgVert(const uiint& child_address) {
        return mPyramidVertChildVert[child_address];
    }
    uiint& getQuadVertProgVert(const uiint& child_address) {
        return mQuadVertChildVert[child_address];
    }
    uiint& getTriangleVertProgVert(const uiint& child_address) {
        return mTriangleVertChildVert[child_address];
    }
    uiint& getBeamVertProgVert(const uiint& child_address) {
        return mBeamVertChildVert[child_address];
    }
    uiint& getVertProgVert(const uiint& ivert, const uiint& elemType);
    uiint& getHexaEdgeProgVert(const uiint& iedge, const uiint& child_address) {
        return mHexaEdgeChildVert[child_address][iedge];
    }
    uiint& getTetraEdgeProgVert(const uiint& iedge, const uiint& child_address) {
        return mTetraEdgeChildVert[child_address][iedge];
    }
    uiint& getPrismEdgeProgVert(const uiint& iedge, const uiint& child_address) {
        return mPrismEdgeChildVert[child_address][iedge];
    }
    uiint& getPyramidEdgeProgVert(const uiint& iedge, const uiint& child_address) {
        return mPyramidEdgeChildVert[child_address][iedge];
    }
    uiint& getQuadEdgeProgVert(const uiint& iedge, const uiint& child_address) {
        return mQuadEdgeChildVert[child_address][iedge];
    }
    uiint& getTriangleEdgeProgVert(const uiint& iedge, const uiint& child_address) {
        return mTriangleEdgeChildVert[child_address][iedge];
    }
    uiint& getBeamEdgeProgVert(const uiint& iedge, const uiint& child_address) {
        return mBeamEdgeChildVert[child_address][iedge];
    }
    uiint& getEdgeProgVert(const uiint& iedge, const uiint& child_address, const uiint& elemType);
    uiint& getHexaFaceProgVert(const uiint& iface, const uiint& child_address) {
        return mHexaFaceChildVert[child_address][iface];
    }
    uiint& getTetraFaceProgVert(const uiint& iface, const uiint& child_address) {
        return mTetraFaceChildVert[child_address][iface];
    }
    uiint& getPrismFaceProgVert(const uiint& iface, const uiint& child_address) {
        return mPrismFaceChildVert[child_address][iface];
    }
    uiint& getPyramidFaceProgVert(const uiint& iface, const uiint& child_address) {
        return mPyramidFaceChildVert[child_address][iface];
    }
    uiint& getQuadFaceProgVert(const uiint& iface, const uiint& child_address) {
        return mQuadFaceChildVert[child_address][iface];
    }
    uiint& getTriangleFaceProgVert(const uiint& iface, const uiint& child_address) {
        return mTriangleFaceChildVert[child_address][iface];
    }
    uiint& getFaceProgVert(const uiint& iface, const uiint& child_address, const uiint& elemType);
    uiint& getHexaVolProgVert(const uiint& child_address) {
        return mHexaVolChildVert[child_address];
    }
    uiint& getTetraVolProgVert(const uiint& child_address) {
        return mTetraVolChildVert[child_address];
    }
    uiint& getPrismVolProgVert(const uiint& child_address) {
        return mPrismVolChildVert[child_address];
    }
    uiint& getPyramidVolProgVert(const uiint& child_address) {
        return mPyramidVolChildVert[child_address];
    }
    uiint& getQuadVolProgVert(const uiint& child_address) {
        return mQuadVolChildVert[child_address];
    }
    uiint& getTriangleVolProgVert(const uiint& child_address) {
        return mTriangleVolChildVert[child_address];
    }
    uiint& getBeamVolProgVert(const uiint& child_address) {
        return mBeamVolChildVert[child_address];
    }
    uiint& getVolProgVert(const uiint& child_address, const uiint& elemType);
};
#endif	/* _PROGELEMENTTREE_H */
}
