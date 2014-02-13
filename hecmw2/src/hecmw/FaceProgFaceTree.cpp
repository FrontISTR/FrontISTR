/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FaceProgFaceTree.cpp
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
#include "FaceProgFaceTree.h"
using namespace pmw;
CFaceProgFace::CFaceProgFace()
{
    uiint nNumOfFace, nNumOfVert;
    uiint iface, ivert;
    CFaceTree *pFaceTree = CFaceTree::Instance();
    nNumOfFace = 6;
    nNumOfVert = 4;
    for(iface=0; iface < nNumOfFace; iface++) {
        for(ivert=0; ivert < nNumOfVert; ivert++) {
            mnHexaProgFace[iface][ivert] = iface;
        };
    };
    nNumOfFace = 4;
    for(iface=0; iface < nNumOfFace; iface++) {
        nNumOfVert = pFaceTree->getPrismFaceNumOfVert(iface);
        uiint* nLocalNode = pFaceTree->getLocalNodeTetraFace(iface);
        for(ivert=0; ivert < nNumOfVert; ivert++) {
            uiint nVertNum = nLocalNode[ivert];
            switch(nVertNum) {
            case(0):
                if(iface==0) mnTetraProgFace[iface][ivert]=0;
                if(iface==1) mnTetraProgFace[iface][ivert]=2;
                if(iface==3) mnTetraProgFace[iface][ivert]=4;
                break;
            case(1):
                if(iface==0) mnTetraProgFace[iface][ivert]=0;
                if(iface==1) mnTetraProgFace[iface][ivert]=2;
                if(iface==2) mnTetraProgFace[iface][ivert]=5;
                break;
            case(2):
                if(iface==0) mnTetraProgFace[iface][ivert]=0;
                if(iface==2) mnTetraProgFace[iface][ivert]=3;
                if(iface==3) mnTetraProgFace[iface][ivert]=4;
                break;
            case(3):
                if(iface==1) mnTetraProgFace[iface][ivert]=2;
                if(iface==2) mnTetraProgFace[iface][ivert]=1;
                if(iface==3) mnTetraProgFace[iface][ivert]=4;
                break;
            }
        };
    };
    nNumOfFace = 5;
    for(iface=0; iface < nNumOfFace; iface++) {
        nNumOfVert = pFaceTree->getPrismFaceNumOfVert(iface);
        uiint* nLocalNode = pFaceTree->getLocalNodePrismFace(iface);
        for(ivert=0; ivert < nNumOfVert; ivert++) {
            uiint nVertNum = nLocalNode[ivert];
            if(iface==0 || iface==1) {
                switch(nVertNum) {
                case(0):
                case(1):
                case(2):
                    mnPrismProgFace[iface][ivert]=0;
                    break;
                case(3):
                case(4):
                case(5):
                    mnPrismProgFace[iface][ivert]=1;
                    break;
                }
            } else {
                switch(nVertNum) {
                case(0):
                    if(iface==2) mnPrismProgFace[iface][ivert]=2;
                    if(iface==4) mnPrismProgFace[iface][ivert]=4;
                    break;
                case(1):
                    if(iface==2) mnPrismProgFace[iface][ivert]=2;
                    if(iface==3) mnPrismProgFace[iface][ivert]=5;
                    break;
                case(2):
                    if(iface==3) mnPrismProgFace[iface][ivert]=3;
                    if(iface==4) mnPrismProgFace[iface][ivert]=4;
                    break;
                case(3):
                    if(iface==2) mnPrismProgFace[iface][ivert]=2;
                    if(iface==4) mnPrismProgFace[iface][ivert]=4;
                    break;
                case(4):
                    if(iface==2) mnPrismProgFace[iface][ivert]=2;
                    if(iface==3) mnPrismProgFace[iface][ivert]=5;
                    break;
                case(5):
                    if(iface==3) mnPrismProgFace[iface][ivert]=3;
                    if(iface==4) mnPrismProgFace[iface][ivert]=4;
                    break;
                }
            }
        };
    };
    nNumOfFace = 1;
    nNumOfVert = 4;
    for(iface=0; iface < nNumOfFace; iface++) {
        for(ivert=0; ivert < nNumOfVert; ivert++) {
            mnQuadProgFace[iface][ivert] = iface;
        };
    };
    nNumOfFace = 1;
    nNumOfVert = 3;
    for(iface=0; iface < nNumOfFace; iface++) {
        for(ivert=0; ivert < nNumOfVert; ivert++) {
            mnTriangleProgFace[iface][ivert] = iface;
        };
    };
}
CFaceProgFace::~CFaceProgFace()
{
    ;
}
uiint& CFaceProgFace::getProgFace_Hexa(const uiint& nFace, const uiint& nFaceVert)
{
    return mnHexaProgFace[nFace][nFaceVert];
}
uiint& CFaceProgFace::getProgFace_Tetra(const uiint& nFace, const uiint& nFaceVert)
{
    return mnTetraProgFace[nFace][nFaceVert];
}
uiint& CFaceProgFace::getProgFace_Prism(const uiint& nFace, const uiint& nFaceVert)
{
    return mnPrismProgFace[nFace][nFaceVert];
}
uiint& CFaceProgFace::getProgFace_Quad(const uiint& nFace, const uiint& nFaceVert)
{
    return mnQuadProgFace[nFace][nFaceVert];
}
uiint& CFaceProgFace::getProgFace_Triangle(const uiint& nFace, const uiint& nFaceVert)
{
    return mnTriangleProgFace[nFace][nFaceVert];
}
