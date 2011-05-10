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
    uint mnHexaProgFace[6][4];//親の面番号-親の面構成頂点番号 =＞ 子の面番号
    uint mnTetraProgFace[4][3];
    uint mnPrismProgFace[5][4];

    uint mnQuadProgFace[1][4];
    uint mnTriangleProgFace[1][3];
    
public:
    uint& getProgFace_Hexa(const uint& nFace, const uint& nFaceVert);
    uint& getProgFace_Tetra(const uint& nFace, const uint& nFaceVert);
    uint& getProgFace_Prism(const uint& nFace, const uint& nFaceVert);
    uint& getProgFace_Quad(const uint& nFace, const uint& nFaceVert);
    uint& getProgFace_Triangle(const uint& nFace, const uint& nFaceVert);
};
#endif	/* FACEPROGFACETREE_H */
}



