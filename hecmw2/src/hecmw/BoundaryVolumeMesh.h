/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryVolumeMesh.h
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
#include "BoundaryVolume.h"
#include "BoundaryMesh.h"
#include "ElementType.h"
#include "ShapeHexa.h"
#include "ShapeTetra.h"
#include "ShapePrism.h"
namespace pmw
{
#ifndef _BOUNDARYVOLUMEMESH_H
#define	_BOUNDARYVOLUMEMESH_H
class CBoundaryVolumeMesh:public CBoundaryMesh
{
public:
    CBoundaryVolumeMesh();
    virtual ~CBoundaryVolumeMesh();
protected:
    vector<CBoundaryVolume*> mvBVolume;
    map<uiint, uiint, less<uiint> > mmBVolumeID2Index;
    vector<CBoundaryNode*> mvBEdgeBNode;
    vector<CBoundaryNode*> mvBFaceBNode;
    vector<CBoundaryNode*> mvBVolBNode;
public:
    void resizeVolume(const uiint& res_size);
    uiint getNumOfVolume() {
        return mvBVolume.size();
    }
    void setBVolume(const uiint& index, CBoundaryVolume *pBVolume);
    void addBVolume(CBoundaryVolume *pBVolume);
    CBoundaryVolume* getBVolumeIX(const uiint& index) {
        return mvBVolume[index];
    }
    CBoundaryVolume* getBVolumeID(const uiint& id) {
        uiint index= mmBVolumeID2Index[id];
        return mvBVolume[index];
    }
    uiint& getBVolumeIndex(const uiint& id) {
        mmBVolumeID2Index[id];
    }
protected:
    vvuint mvAggregateVol;
public:
    void resizeAggVol();
    void setupAggVol();
    void setAggVol(const uiint& ibnode, const uiint& nVolumeID) {
        mvAggregateVol[ibnode].push_back(nVolumeID);
    }
    void GeneEdgeBNode();
    void GeneFaceBNode();
    void GeneVolBNode();
    void refine(CBoundaryVolumeMesh *pProgBVolMesh);
    void deleteProgData();
protected:
    virtual void distNeumannValue();
    virtual void distDirichletValue();
};
#endif	/* _BOUNDARYVOLUMEMESH_H */
}
