/* 
 * File:   Quad2.h
 * Author: ktakeda
 *
 * Created on 2010/11/19, 14:59
 */
#include "Quad.h"

namespace pmw{
#ifndef QUAD2_H
#define	QUAD2_H
class CQuad2:public CQuad{
public:
    CQuad2();
    virtual ~CQuad2();

private:
    static uint mnElemType;
    static uint mnElemOrder;
    static uint mNumOfFace;
    static uint mNumOfEdge;
    static uint mNumOfNode;
    static uint mNumOfVert;

public:
    virtual void initialize();

public:
    // Property
    virtual const uint& getType(){ return mnElemType;}
    virtual const uint& getOrder(){ return mnElemOrder;}
    virtual const uint& getNumOfFace(){ return mNumOfFace;}
    virtual const uint& getNumOfEdge(){ return mNumOfEdge;}
    virtual const uint& getNumOfNode(){ return mNumOfNode;}
    virtual const uint& getNumOfVert(){ return mNumOfVert;}

    //
    // 2次要素において、辺NodeをmvNodeに移し替える && 1次要素では何もしない.
    //
    virtual void replaseEdgeNode();

//    // Refine後
//    // 1. 辺-面 Element*配列を解放
//    // 2. 辺-面 Node*   配列を解放 (2次要素は辺ノードを残す)
//    // --
//    virtual void deleteProgData();
};
#endif	/* QUAD2_H */
}


