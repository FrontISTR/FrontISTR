//
//  Vertex.h
//
//                          2009.09.24
//                          2008.12.18
//                          k.Takeda
#ifndef VERTEX_HH_E70B6426_2710_4010_95B0_6219107F5EDF
#define VERTEX_HH_E70B6426_2710_4010_95B0_6219107F5EDF

#include "CommonStd.h"
#include "TypeDef.h"

#include <map>

namespace pmw{
class CVertex{
public:
    CVertex(void);
    virtual ~CVertex(void);

protected:
    uint     mID;//uint mIndex;
    vdouble mvCoord;
    
    vuint mvAggElementID;   //頂点を共有している要素番号(ID)配列(AggElementのID)
    map<uint, uint, less<uint> > mmNeibElemVertNum;//接続している隣の"要素ID -> 頂点番号"

public:
    // Node ID == Node Index(in Mesh)
    void setID(const uint& id){ mID = id;}
    uint& getID(){ return mID;}
    // 
    //void setIndex(const uint& index){ mIndex=index;}
    //uint& getID(){ return mIndex;}


    // 座標
    void     setCoord(const vdouble& coord){ mvCoord = coord;}
    vdouble& getCoord(){ return mvCoord;}

    double& getX(){ return mvCoord[0];}
    double& getY(){ return mvCoord[1];}
    double& getZ(){ return mvCoord[2];}

    // for AggregateElement
    void  setAggElemID(const uint& id){ mvAggElementID.push_back(id);}
    uint  getNumOfAggElem(){ return  mvAggElementID.size();} // Vertexが所属する要素の個数
    uint& getAggElemID(const uint& i){ return mvAggElementID[i];}
    void  clearAggElemID(){ mvAggElementID.clear();}// 2段目以降のprolongationのために必要.

    // for AggregateElement 追加(隣の要素の頂点番号) <= CommMeshのIndex管理で使用
    void setNeibElemVert(const uint& elemID,const uint& localNum){ mmNeibElemVertNum[elemID]=localNum;}
    uint getNumOfNeibElemVert(){ return mmNeibElemVertNum.size();}
    uint& getNeibElemIDVert(const uint& elemID){ return mmNeibElemVertNum[elemID];}
    void clearNeibElemVert(){ mmNeibElemVertNum.clear();}// 2段目以降のprolongationのために必要.

};
}
#endif
