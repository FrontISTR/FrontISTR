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
    uint mID;

    vuint mvAggElementID; //頂点を共有している要素番号(ID)配列
    map<uint, uint, less<uint> > mmNeibElemVertNum;//接続している隣の"要素ID -> 頂点番号"

public:
    // Node ID
    void setID(const uint& id){ mID = id;}
    uint& getID(){ return mID;}
    
    // AggregateElement
    // ----
    void  setAggElemID(const uint& id){ mvAggElementID.push_back(id);}
    uint  getNumOfAggElem(){ return  mvAggElementID.size();} // Vertexが所属する要素の個数
    uint& getAggElemID(const uint& i){ return mvAggElementID[i];}
    void  clearAggElemID(){ mvAggElementID.clear();}// 2段目以降のprolongationのために必要.

    // AggregateElement 追加(隣の要素の頂点番号) <= CommMeshのIndex管理で使用
    // ----
    void setNeibElemVert(const uint& elemID,const uint& localNum){ mmNeibElemVertNum[elemID]=localNum;}
    uint getNumOfNeibElemVert(){ return mmNeibElemVertNum.size();}
    uint& getNeibElemIDVert(const uint& elemID){ return mmNeibElemVertNum[elemID];}
    void clearNeibElemVert(){ mmNeibElemVertNum.clear();}// 2段目以降のprolongationのために必要.
};
#endif	/* _BVERTEX_H */
}

