/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BndVertex.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "TypeDef.h"
#include <map>
namespace pmw
{
#ifndef _BndVERTEX_H
#define	_BndVERTEX_H
class CBndVertex
{
public:
    CBndVertex();
    virtual ~CBndVertex();
protected:
    uiint mID;
    vuint mvAggElementID;
    map<uiint, uiint, less<uiint> > mmNeibElemVertNum;
public:
    void setID(const uiint& id){ mID = id;}
    uiint& getID(){ return mID;}
    
    void  setAggElemID(const uiint& id){ mvAggElementID.push_back(id);}
    uiint  getNumOfAggElem(){ return  mvAggElementID.size();} 
    uiint& getAggElemID(const uiint& i){ return mvAggElementID[i];}
    void  clearAggElemID(){ mvAggElementID.clear();}

    void setNeibElemVert(const uiint& elemID,const uiint& localNum){ mmNeibElemVertNum[elemID]=localNum;}
    uiint getNumOfNeibElemVert(){ return mmNeibElemVertNum.size();}
    uiint& getNeibElemIDVert(const uiint& elemID){ return mmNeibElemVertNum[elemID];}
    void clearNeibElemVert(){ mmNeibElemVertNum.clear();}
    void deleteAggregate();
};
#endif	/* _BVERTEX_H */
}
