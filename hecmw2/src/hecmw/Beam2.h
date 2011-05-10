/* 
 * File:   Beam2.h
 * Author: ktakeda
 *
 * Created on 2010/11/19, 15:00
 */
#include "Beam.h"


namespace pmw{
#ifndef BEAM2_H
#define	BEAM2_H
class CBeam2:public CBeam{
public:
    CBeam2();
    virtual ~CBeam2();

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
#endif	/* BEAM2_H */
}

