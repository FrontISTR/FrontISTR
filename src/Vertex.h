//
//  Vertex.h
//
//                          2008.12.18
//                          2008.12.18
//                          k.Takeda
#ifndef VERTEX_HH_E70B6426_2710_4010_95B0_6219107F5EDF
#define VERTEX_HH_E70B6426_2710_4010_95B0_6219107F5EDF

#include "CommonStd.h"
#include "TypeDef.h"

namespace pmw{
class CVertex{
public:
    CVertex(void);
    virtual ~CVertex(void);

protected:
    int     mID;//uint mIndex;
    vdouble mvCoord;

    vuint mvElementIndex;//頂点が屬している要素番号(Index)

public:
    // Node ID == Node Index(in Mesh)
    void setID(const int& id){ mID = id;}
    int& getID(){ return mID;}
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
    void  setElemIndex(const uint& index){ mvElementIndex.push_back(index);}
    uint  getNumOfAggElem(){ return  mvElementIndex.size();} // Vertexが所属する要素の個数
    uint& getElemIndex(const uint& i){ return mvElementIndex[i];}

    void  clearElemIndex(){ mvElementIndex.clear();}//2段 以降のprolongationのために必要.

};
}
#endif
