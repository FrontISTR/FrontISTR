/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BNodeMeshGrp.h
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
#include "BoundaryGroup.h"
#include "BoundaryNodeMesh.h"
namespace pmw
{
#ifndef _BNODEMESH_GRP_H
#define	_BNODEMESH_GRP_H
class CBNodeMeshGrp
{
public:
    CBNodeMeshGrp();
    virtual ~CBNodeMeshGrp();
protected:
    BoundaryGroup<CBoundaryNodeMesh*> mGrpBndNodeMesh;
public:
    void reserveBndNodeMesh(const uiint& res_size) {
        mGrpBndNodeMesh.reserve(res_size);
    }
    void setBndNodeMesh(CBoundaryNodeMesh *pBNodeMesh) {
        mGrpBndNodeMesh.push(pBNodeMesh);
    }
    CBoundaryNodeMesh* getBndNodeMeshIX(const uiint& index) {
        return mGrpBndNodeMesh.get_withIndex(index);
    }
    CBoundaryNodeMesh* getBndNodeMeshID(const uiint& id) {
        return mGrpBndNodeMesh.get_withID(id);
    }
    uiint getNumOfBoundaryNodeMesh() {
        return mGrpBndNodeMesh.NumOfBoundary();
    }

    void clear();
};
#endif	/* _BNODEMESH_GRP_H */
}
