/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BoundaryNode.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _BOUNDARYNODE_H_37946911_224e_49ea_90c9_472937aab367
#define	_BOUNDARYNODE_H_37946911_224e_49ea_90c9_472937aab367
#include "Vertex.h"
#include "GeneralBoundary.h"
#include "Logger.h"
namespace pmw{
class CBoundaryNode:public CVertex, public CGeneralBoundary{
public:
    CBoundaryNode();
    virtual ~CBoundaryNode();
private:
public:
    virtual void initialize(const uint& bndType,const uint& dof);
};
}
#endif	/* _BOUNDARYNODE_H */
