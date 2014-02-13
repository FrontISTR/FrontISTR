/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Quad2.h
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
#include "Quad.h"
namespace pmw
{
#ifndef QUAD2_H
#define	QUAD2_H
class CQuad2:public CQuad
{
public:
    CQuad2();
    virtual ~CQuad2();
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
#endif	/* QUAD2_H */
}
