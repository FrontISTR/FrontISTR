/* 
 * File:   BndVertex.h
 * Author: ktakeda
 *
 * 座標を持たないVertex
 *
 * Vertexの親
 *
 * Created on 2010/04/28, 12:09
 */
#include "TypeDef.h"
#include <map>

namespace pmw{
#ifndef _BndVERTEX_H
#define	_BndVERTEX_H
class CBndVertex{
public:
    CBndVertex();
    virtual ~CBndVertex();

protected:
    uiint mID;

    vuint mvAggElementID; //頂点を共有している要素番号(ID)配列
    map<uiint, uiint, less<uiint> > mmNeibElemVertNum;//接続している隣の"要素ID -> 頂点番号"

public:
    // Node ID
    void setID(const uiint& id){ mID = id;}
    uiint& getID(){ return mID;}
    
    // AggregateElement
    // ----
    void  setAggElemID(const uiint& id){ mvAggElementID.push_back(id);}
    uiint  getNumOfAggElem(){ return  mvAggElementID.size();} // Vertexが所属する要素の個数
    uiint& getAggElemID(const uiint& i){ return mvAggElementID[i];}
    void  clearAggElemID(){ mvAggElementID.clear();}// 2段目以降のprolongationのために必要.

    // AggregateElement 追加(隣の要素の頂点番号) <= CommMeshのIndex管理で使用
    // ----
    void setNeibElemVert(const uiint& elemID,const uiint& localNum){ mmNeibElemVertNum[elemID]=localNum;}
    uiint getNumOfNeibElemVert(){ return mmNeibElemVertNum.size();}
    uiint& getNeibElemIDVert(const uiint& elemID){ return mmNeibElemVertNum[elemID];}
    void clearNeibElemVert(){ mmNeibElemVertNum.clear();}// 2段目以降のprolongationのために必要.
    
    
    // AggregateElementの解放
    //
    // |全てのMesh処理が終わった後で呼び出す| = CMW::FinalizeRefine()の処理
    //
    void deleteAggregate();
};
#endif	/* _BVERTEX_H */
}

