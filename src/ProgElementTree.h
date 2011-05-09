/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ProgElementTree.h
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
#include "ElementType.h"
#include "Logger.h"
namespace pmw{
#ifndef _PROGELEMENTTREE_H
#define	_PROGELEMENTTREE_H
class CProgElementTree{
private:
    CProgElementTree();
public:
    static CProgElementTree* Instance(){
        static CProgElementTree progElemTree;
        return &progElemTree;
    }
    virtual ~CProgElementTree();
protected:
    uint mInvalidNum;
    uint mProgVert;
    uint mHexaVertChildVert[8];
    uint mTetraVertChildVert[4];
    uint mPrismVertChildVert[6];
    uint mPyramidVertChildVert[5];
    uint mQuadVertChildVert[4];
    uint mTriangleVertChildVert[3];
    uint mBeamVertChildVert[2];
    uint mHexaEdgeChildVert[8][12];
    uint mTetraEdgeChildVert[4][6];
    uint mPrismEdgeChildVert[6][9];
    uint mPyramidEdgeChildVert[8][8];
    uint mQuadEdgeChildVert[4][4];
    uint mTriangleEdgeChildVert[3][3];
    uint mBeamEdgeChildVert[2][1];
    uint mHexaFaceChildVert[8][6];
    uint mTetraFaceChildVert[4][4];
    uint mPrismFaceChildVert[6][5];
    uint mPyramidFaceChildVert[8][5];
    uint mQuadFaceChildVert[4][1];
    uint mTriangleFaceChildVert[3][1];
    uint mHexaVolChildVert[8];
    uint mTetraVolChildVert[4];
    uint mPrismVolChildVert[6];
    uint mPyramidVolChildVert[8];
    uint mQuadVolChildVert[4];
    uint mTriangleVolChildVert[3];
    uint mBeamVolChildVert[2];
public:
    uint& getInvalidNum(){ return mInvalidNum;}
    uint& getHexaVertProgVert(const uint& child_address){ return mHexaVertChildVert[child_address];}
    uint& getTetraVertProgVert(const uint& child_address){ return mTetraVertChildVert[child_address];}
    uint& getPrismVertProgVert(const uint& child_address){ return mPrismVertChildVert[child_address];}
    uint& getPyramidVertProgVert(const uint& child_address){ return mPyramidVertChildVert[child_address];}
    uint& getQuadVertProgVert(const uint& child_address){ return mQuadVertChildVert[child_address];}
    uint& getTriangleVertProgVert(const uint& child_address){ return mTriangleVertChildVert[child_address];}
    uint& getBeamVertProgVert(const uint& child_address){ return mBeamVertChildVert[child_address];}
    uint& getVertProgVert(const uint& ivert, const uint& elemType);
    uint& getHexaEdgeProgVert(const uint& iedge, const uint& child_address){ return mHexaEdgeChildVert[child_address][iedge];}
    uint& getTetraEdgeProgVert(const uint& iedge, const uint& child_address){ return mTetraEdgeChildVert[child_address][iedge];}
    uint& getPrismEdgeProgVert(const uint& iedge, const uint& child_address){ return mPrismEdgeChildVert[child_address][iedge];}
    uint& getPyramidEdgeProgVert(const uint& iedge, const uint& child_address){ return mPyramidEdgeChildVert[child_address][iedge];}
    uint& getQuadEdgeProgVert(const uint& iedge, const uint& child_address){ return mQuadEdgeChildVert[child_address][iedge];}
    uint& getTriangleEdgeProgVert(const uint& iedge, const uint& child_address){ return mTriangleEdgeChildVert[child_address][iedge];}
    uint& getBeamEdgeProgVert(const uint& iedge, const uint& child_address){ return mBeamEdgeChildVert[child_address][iedge];}
    uint& getEdgeProgVert(const uint& iedge, const uint& child_address, const uint& elemType);
    uint& getHexaFaceProgVert(const uint& iface, const uint& child_address){ return mHexaFaceChildVert[child_address][iface];}
    uint& getTetraFaceProgVert(const uint& iface, const uint& child_address){ return mTetraFaceChildVert[child_address][iface];}
    uint& getPrismFaceProgVert(const uint& iface, const uint& child_address){ return mPrismFaceChildVert[child_address][iface];}
    uint& getPyramidFaceProgVert(const uint& iface, const uint& child_address){ return mPyramidFaceChildVert[child_address][iface];}
    uint& getQuadFaceProgVert(const uint& iface, const uint& child_address){ return mQuadFaceChildVert[child_address][iface];}
    uint& getTriangleFaceProgVert(const uint& iface, const uint& child_address){ return mTriangleFaceChildVert[child_address][iface];}
    uint& getFaceProgVert(const uint& iface, const uint& child_address, const uint& elemType);
    uint& getHexaVolProgVert(const uint& child_address){ return mHexaVolChildVert[child_address];}
    uint& getTetraVolProgVert(const uint& child_address){ return mTetraVolChildVert[child_address];}
    uint& getPrismVolProgVert(const uint& child_address){ return mPrismVolChildVert[child_address];}
    uint& getPyramidVolProgVert(const uint& child_address){ return mPyramidVolChildVert[child_address];}
    uint& getQuadVolProgVert(const uint& child_address){ return mQuadVolChildVert[child_address];}
    uint& getTriangleVolProgVert(const uint& child_address){ return mTriangleVolChildVert[child_address];}
    uint& getBeamVolProgVert(const uint& child_address){ return mBeamVolChildVert[child_address];}
    uint& getVolProgVert(const uint& child_address, const uint& elemType);
};
#endif	/* _PROGELEMENTTREE_H */
}
