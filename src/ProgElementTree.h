/* 
 * File:   ProgElementTree.h
 * Author: ktakeda
 *
 * prolongation時の親Elementの辺,面の点と,子Element(progElem)の頂点の関係
 *
 * Created on 2009/09/03, 17:23
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

    uint mProgVert;//汎用メソッドの返り値用メンバー

    //{親の頂点} => 子(progElem)の頂点
    uint mHexaVertChildVert[8];
    uint mTetraVertChildVert[4];
    uint mPrismVertChildVert[6];
    uint mPyramidVertChildVert[5];
    uint mQuadVertChildVert[4];
    uint mTriangleVertChildVert[3];
    uint mBeamVertChildVert[2];

    //{親の辺番号} => 子(progElem)の頂点番号:子(progElem)のAddress(親の頂点番号)別
    uint mHexaEdgeChildVert[8][12];
    uint mTetraEdgeChildVert[4][6];
    uint mPrismEdgeChildVert[6][9];
    uint mPyramidEdgeChildVert[8][8];//5番めの頂点には4個のPyramidがぶら下がる
    uint mQuadEdgeChildVert[4][4];
    uint mTriangleEdgeChildVert[3][3];
    uint mBeamEdgeChildVert[2][1];

    //{親の面番号} => 子(progElem)の頂点番号:子(progElem)のAddress(親の頂点番号)別
    uint mHexaFaceChildVert[8][6];
    uint mTetraFaceChildVert[4][4];
    uint mPrismFaceChildVert[6][5];
    uint mPyramidFaceChildVert[8][5];//5番めの頂点には4個のPyramidがぶら下がる
    uint mQuadFaceChildVert[4][1];
    uint mTriangleFaceChildVert[3][1];
    
    //{要素中心} => 子(progElem)の頂点番号:子(progElem)のAddress(親の頂点番号)別
    uint mHexaVolChildVert[8];
    uint mTetraVolChildVert[4];
    uint mPrismVolChildVert[6];
    uint mPyramidVolChildVert[8];//5番めの頂点には4個のPyramidがぶら下がる
    uint mQuadVolChildVert[4];
    uint mTriangleVolChildVert[3];
    uint mBeamVolChildVert[2];

public:
    // 子供の頂点番号が無い場合の数値.
    uint& getInvalidNum(){ return mInvalidNum;}

    // 子(progElem)のアドレス(親の頂点番号)から,子(progElem)の頂点番号を提供
    // --
    uint& getHexaVertProgVert(const uint& child_address){ return mHexaVertChildVert[child_address];}
    uint& getTetraVertProgVert(const uint& child_address){ return mTetraVertChildVert[child_address];}
    uint& getPrismVertProgVert(const uint& child_address){ return mPrismVertChildVert[child_address];}
    uint& getPyramidVertProgVert(const uint& child_address){ return mPyramidVertChildVert[child_address];}
    uint& getQuadVertProgVert(const uint& child_address){ return mQuadVertChildVert[child_address];}
    uint& getTriangleVertProgVert(const uint& child_address){ return mTriangleVertChildVert[child_address];}
    uint& getBeamVertProgVert(const uint& child_address){ return mBeamVertChildVert[child_address];}

    uint& getVertProgVert(const uint& ivert, const uint& elemType);


    // 子(progElem)のアドレス(親の頂点番号)と,親(Elem)の辺番号から,子供(progElem)の頂点番号を提供
    // --
    // iedge:親の辺番号, child_address:子供が居る場所(親の頂点番号)
    uint& getHexaEdgeProgVert(const uint& iedge, const uint& child_address){ return mHexaEdgeChildVert[child_address][iedge];}
    uint& getTetraEdgeProgVert(const uint& iedge, const uint& child_address){ return mTetraEdgeChildVert[child_address][iedge];}
    uint& getPrismEdgeProgVert(const uint& iedge, const uint& child_address){ return mPrismEdgeChildVert[child_address][iedge];}
    uint& getPyramidEdgeProgVert(const uint& iedge, const uint& child_address){ return mPyramidEdgeChildVert[child_address][iedge];}
    uint& getQuadEdgeProgVert(const uint& iedge, const uint& child_address){ return mQuadEdgeChildVert[child_address][iedge];}
    uint& getTriangleEdgeProgVert(const uint& iedge, const uint& child_address){ return mTriangleEdgeChildVert[child_address][iedge];}
    uint& getBeamEdgeProgVert(const uint& iedge, const uint& child_address){ return mBeamEdgeChildVert[child_address][iedge];}

    uint& getEdgeProgVert(const uint& iedge, const uint& child_address, const uint& elemType);


    // 子(progElem)のアドレス(親の頂点番号)と,親(Elem)の面番号から,子供(progElem)の頂点番号を提供
    // --
    // iface:親の面番号, child_address:子供が居る場所(親の頂点番号)
    uint& getHexaFaceProgVert(const uint& iface, const uint& child_address){ return mHexaFaceChildVert[child_address][iface];}
    uint& getTetraFaceProgVert(const uint& iface, const uint& child_address){ return mTetraFaceChildVert[child_address][iface];}
    uint& getPrismFaceProgVert(const uint& iface, const uint& child_address){ return mPrismFaceChildVert[child_address][iface];}
    uint& getPyramidFaceProgVert(const uint& iface, const uint& child_address){ return mPyramidFaceChildVert[child_address][iface];}
    uint& getQuadFaceProgVert(const uint& iface, const uint& child_address){ return mQuadFaceChildVert[child_address][iface];}
    uint& getTriangleFaceProgVert(const uint& iface, const uint& child_address){ return mTriangleFaceChildVert[child_address][iface];}

    uint& getFaceProgVert(const uint& iface, const uint& child_address, const uint& elemType);



    // 子(progElem)のアドレス(親の頂点番号)から,親の要素中心 => 子供(progElem)の頂点番号を提供
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

