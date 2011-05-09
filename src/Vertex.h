/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Vertex.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
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
    uint    mID;
    vdouble mvCoord;
    vuint mvAggElementID; 
    map<uint, uint, less<uint> > mmNeibElemVertNum;
public:
    void setID(const uint& id){ mID = id;}
    uint& getID(){ return mID;}
    void     setCoord(const vdouble& coord){ mvCoord = coord;}
    vdouble& getCoord(){ return mvCoord;}
    double& getX(){ return mvCoord[0];}
    double& getY(){ return mvCoord[1];}
    double& getZ(){ return mvCoord[2];}
    void  setAggElemID(const uint& id){ mvAggElementID.push_back(id);}
    uint  getNumOfAggElem(){ return  mvAggElementID.size();} 
    uint& getAggElemID(const uint& i){ return mvAggElementID[i];}
    void  clearAggElemID(){ mvAggElementID.clear();}
    void setNeibElemVert(const uint& elemID,const uint& localNum){ mmNeibElemVertNum[elemID]=localNum;}
    uint getNumOfNeibElemVert(){ return mmNeibElemVertNum.size();}
    uint& getNeibElemIDVert(const uint& elemID){ return mmNeibElemVertNum[elemID];}
    void clearNeibElemVert(){ mmNeibElemVertNum.clear();}
};
}
#endif
