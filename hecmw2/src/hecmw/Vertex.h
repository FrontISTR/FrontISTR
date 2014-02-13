/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Vertex.h
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
#include "BndVertex.h"
namespace pmw
{
#ifndef _VERTEX_HH_E70B6426_
#define _VERTEX_HH_E70B6426_
class CVertex:public CBndVertex
{
public:
    CVertex(void);
    virtual ~CVertex(void);
protected:
    vdouble mvCoord;
public:
    void     setCoord(const vdouble& coord) {
        mvCoord = coord;
    }
    vdouble& getCoord() {
        return mvCoord;
    }
    double& getX() {
        return mvCoord[0];
    }
    double& getY() {
        return mvCoord[1];
    }
    double& getZ() {
        return mvCoord[2];
    }
};
#endif /* _VERTEX_HH_E70B6426_  */
}
