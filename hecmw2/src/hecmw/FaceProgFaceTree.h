/* 
 * File:   FaceProgFaceTree.h
 * Author: ktakeda
 *
 * Created on 2010/10/15, 21:11
 */
#include "CommonStd.h"
#include "TypeDef.h"
#include "FaceTree.h"

namespace pmw{
#ifndef FACEPROGFACETREE_H
#define	FACEPROGFACETREE_H
class CFaceProgFace{
private:
    CFaceProgFace();
public:
    static CFaceProgFace* Instance(){
        static CFaceProgFace moProgFace;
        return &moProgFace;
    }
    ~CFaceProgFace();

private:
    uiint mnHexaProgFace[6][4];//親の面番号-親の面構成頂点番号 =＞ 子の面番号
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



