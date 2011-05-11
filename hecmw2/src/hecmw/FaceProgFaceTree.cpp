//
// FaceProgFaceTree.cpp
//
//              2010.10.15
//              k.Takeda
#include "FaceProgFaceTree.h"
using namespace pmw;


CFaceProgFace::CFaceProgFace()
{
    uiint nNumOfFace, nNumOfVert;
    uiint iface, ivert;

    CFaceTree *pFaceTree = CFaceTree::Instance();

    /*
     SkinFace::refine()
    */
    // ------------------------------------
    // [親の面番号][親の頂点番号] =＞ 子の面番号
    // ------------------------------------
    // 1.Hexa 分割しても面番号は変わらない.
    //
    nNumOfFace = 6; nNumOfVert = 4;
    for(iface=0; iface < nNumOfFace; iface++){
        for(ivert=0; ivert < nNumOfVert; ivert++){
            mnHexaProgFace[iface][ivert] = iface;
        };
    };
    // 
    // 2. Tetra (分割した要素は6面体)
    //
    nNumOfFace = 4;
    for(iface=0; iface < nNumOfFace; iface++){

        nNumOfVert = pFaceTree->getPrismFaceNumOfVert(iface);      //面を構成する節点数
        uiint* nLocalNode = pFaceTree->getLocalNodeTetraFace(iface);//面を構成する頂点の”要素構成 頂点番号”

        for(ivert=0; ivert < nNumOfVert; ivert++){

            uiint nVertNum = nLocalNode[ivert];//要素構成の頂点番号

            switch(nVertNum){
                case(0):
                    if(iface==0) mnTetraProgFace[iface][ivert]=0;//(0); // <- 親の面番号=0 頂点番号=nVertNum に位置する子要素の"親の面=0に接する"面番号
                    if(iface==1) mnTetraProgFace[iface][ivert]=2;//(2);
                    if(iface==3) mnTetraProgFace[iface][ivert]=4;//(4);
                    break;
                case(1):
                    if(iface==0) mnTetraProgFace[iface][ivert]=0;//(0);
                    if(iface==1) mnTetraProgFace[iface][ivert]=2;//(2);
                    if(iface==2) mnTetraProgFace[iface][ivert]=5;//(5);
                    break;
                case(2):
                    if(iface==0) mnTetraProgFace[iface][ivert]=0;//(0);
                    if(iface==2) mnTetraProgFace[iface][ivert]=3;//(3);
                    if(iface==3) mnTetraProgFace[iface][ivert]=4;//(4);
                    break;
                case(3):
                    if(iface==1) mnTetraProgFace[iface][ivert]=2;//(2);
                    if(iface==2) mnTetraProgFace[iface][ivert]=1;//(1);
                    if(iface==3) mnTetraProgFace[iface][ivert]=4;//(4);
                    break;
            }//switch(頂点番号)
        };
    };
    // 3. Prism
    //
    nNumOfFace = 5;
    for(iface=0; iface < nNumOfFace; iface++){

        nNumOfVert = pFaceTree->getPrismFaceNumOfVert(iface);      //面を構成する節点数
        uiint* nLocalNode = pFaceTree->getLocalNodePrismFace(iface);//面を構成する頂点の”要素構成 頂点番号”

        for(ivert=0; ivert < nNumOfVert; ivert++){
            
            uiint nVertNum = nLocalNode[ivert];

            // Prismの三角形の面の場合(面番号==0,1)
            // * 頂点 0,1,2 => 面番号 0
            // * 頂点 3,4,5 => 面番号 1
            if(iface==0 || iface==1){
                switch(nVertNum){
                    case(0):case(1):case(2):
                        mnPrismProgFace[iface][ivert]=0;//(0);
                        break;
                    case(3):case(4):case(5):
                        mnPrismProgFace[iface][ivert]=1;//(1);
                        break;
                }//switch(頂点番号)
            }else{
                switch(nVertNum){
                    case(0):
                        if(iface==2) mnPrismProgFace[iface][ivert]=2;//(2); // <- 親の面番号=2 頂点番号=nVertNum に位置する子要素の, "親の面=2に接する"面番号
                        if(iface==4) mnPrismProgFace[iface][ivert]=4;//(4);
                        break;
                    case(1):
                        if(iface==2) mnPrismProgFace[iface][ivert]=2;//(2);
                        if(iface==3) mnPrismProgFace[iface][ivert]=5;//(5);
                        break;
                    case(2):
                        if(iface==3) mnPrismProgFace[iface][ivert]=3;//(3);
                        if(iface==4) mnPrismProgFace[iface][ivert]=4;//(4);
                        break;
                    case(3):
                        if(iface==2) mnPrismProgFace[iface][ivert]=2;//(2);
                        if(iface==4) mnPrismProgFace[iface][ivert]=4;//(4);
                        break;
                    case(4):
                        if(iface==2) mnPrismProgFace[iface][ivert]=2;//(2);
                        if(iface==3) mnPrismProgFace[iface][ivert]=5;//(5);
                        break;
                    case(5):
                        if(iface==3) mnPrismProgFace[iface][ivert]=3;//(3);
                        if(iface==4) mnPrismProgFace[iface][ivert]=4;//(4);
                        break;
                }//switch(頂点番号)
            }//endif
        };
    };
        

    // 4.Quad 分割しても面番号は変わらない.
    //
    nNumOfFace = 1; nNumOfVert = 4;
    for(iface=0; iface < nNumOfFace; iface++){
        for(ivert=0; ivert < nNumOfVert; ivert++){
            mnQuadProgFace[iface][ivert] = iface;
        };
    };
    // 5.Triangle 分割しても面番号は変わらない.
    //
    nNumOfFace = 1; nNumOfVert = 3;
    for(iface=0; iface < nNumOfFace; iface++){
        for(ivert=0; ivert < nNumOfVert; ivert++){
            mnTriangleProgFace[iface][ivert] = iface;
        };
    };
}
CFaceProgFace::~CFaceProgFace()
{
    ;
}

// 親の面番号-親の面構成頂点番号 =＞ 子の面番号
// 
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



