/* 
 * File:   Prism2.h
 * Author: ktakeda
 *
 * Created on 2010/11/19, 14:58
 */
#include "Prism.h"

namespace pmw{
#ifndef PRISM2_H
#define	PRISM2_H
class CPrism2:public CPrism{
public:
    CPrism2();
    virtual ~CPrism2();

private:
    static uiint mnElemType;
    static uiint mnElemOrder;
    static uiint mNumOfFace;
    static uiint mNumOfEdge;
    static uiint mNumOfNode;
    static uiint mNumOfVert;

public:
    virtual void initialize();

public:
    // Property
    virtual const uiint& getType(){ return mnElemType;}
    virtual const uiint& getOrder(){ return mnElemOrder;}
    virtual const uiint& getNumOfFace(){ return mNumOfFace;}
    virtual const uiint& getNumOfEdge(){ return mNumOfEdge;}
    virtual const uiint& getNumOfNode(){ return mNumOfNode;}
    virtual const uiint& getNumOfVert(){ return mNumOfVert;}

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
#endif	/* PRISM2_H */
}

