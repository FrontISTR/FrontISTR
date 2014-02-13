/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Tetra2.h
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
#include "Tetra.h"
namespace pmw
{
#ifndef TETRA2_H
#define	TETRA2_H
class CTetra2:public CTetra
{
public:
    CTetra2();
    virtual ~CTetra2();
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
    virtual const uiint& getType() {
        return mnElemType;
    }
    virtual const uiint& getOrder() {
        return mnElemOrder;
    }
    virtual const uiint& getNumOfFace() {
        return mNumOfFace;
    }
    virtual const uiint& getNumOfEdge() {
        return mNumOfEdge;
    }
    virtual const uiint& getNumOfNode() {
        return mNumOfNode;
    }
    virtual const uiint& getNumOfVert() {
        return mNumOfVert;
    }
    virtual void replaseEdgeNode();
};
#endif	/* TETRA2_H */
}
