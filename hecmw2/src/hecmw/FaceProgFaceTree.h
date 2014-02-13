/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FaceProgFaceTree.h
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
#include "CommonStd.h"
#include "TypeDef.h"
#include "FaceTree.h"
namespace pmw
{
#ifndef FACEPROGFACETREE_H
#define	FACEPROGFACETREE_H
class CFaceProgFace
{
private:
    CFaceProgFace();
public:
    static CFaceProgFace* Instance() {
        static CFaceProgFace moProgFace;
        return &moProgFace;
    }
    ~CFaceProgFace();
private:
    uiint mnHexaProgFace[6][4];
    uiint mnTetraProgFace[4][3];
    uiint mnPrismProgFace[5][4];
    uiint mnQuadProgFace[1][4];
    uiint mnTriangleProgFace[1][3];
public:
    uiint& getProgFace_Hexa(const uiint& nFace, const uiint& nFaceVert);
    uiint& getProgFace_Tetra(const uiint& nFace, const uiint& nFaceVert);
    uiint& getProgFace_Prism(const uiint& nFace, const uiint& nFaceVert);
    uiint& getProgFace_Quad(const uiint& nFace, const uiint& nFaceVert);
    uiint& getProgFace_Triangle(const uiint& nFace, const uiint& nFaceVert);
};
#endif	/* FACEPROGFACETREE_H */
}
