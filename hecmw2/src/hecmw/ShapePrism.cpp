//
//  ShapePrism.cpp
//
//
//
//
//              2010.01.29
//              k.Takeda
#include "ShapePrism.h"
using namespace pmw;


CShapePrism::CShapePrism()
{
//    // MW3 辺番号 -> FrontISTR 2次節点番号
//    mvEdge_2_ISTR[0]=  8;
//    mvEdge_2_ISTR[1]=  7;
//    mvEdge_2_ISTR[2]=  6;
//    mvEdge_2_ISTR[3]= 12;
//    mvEdge_2_ISTR[4]= 13;
//    mvEdge_2_ISTR[5]= 14;
//    mvEdge_2_ISTR[6]= 11;
//    mvEdge_2_ISTR[7]=  9;
//    mvEdge_2_ISTR[8]= 10;
//
//    // FrontISTR -> MW3 辺番号
//    uint i;
//    for(i=0; i< 6; i++){
//        mvISTR_2_Edge[i]= i;
//    };
//    mvISTR_2_Edge[6]=  2;
//    mvISTR_2_Edge[7]=  1;
//    mvISTR_2_Edge[8]=  0;
//    mvISTR_2_Edge[9]=  7;
//    mvISTR_2_Edge[10]= 8;
//    mvISTR_2_Edge[11]= 6;
//    mvISTR_2_Edge[12]= 3;
//    mvISTR_2_Edge[13]= 4;
//    mvISTR_2_Edge[14]= 5;

    //--
    //積分座標(面積座標,gzi): L1,L2,gzi
    //--
    // 2点 積分
    mGzi2[0][0]= 0.333333333333333; mGzi2[0][1]= 0.333333333333333; mGzi2[0][2]= -0.577350269189626;
    mGzi2[1][0]= 0.333333333333333; mGzi2[1][1]= 0.333333333333333; mGzi2[1][2]=  0.577350269189626;
    // 6点 積分
    mGzi6[0][0]= 0.166666666666667; mGzi6[0][1]= 0.166666666666667; mGzi6[0][2]= -0.577350269189626;
    mGzi6[1][0]= 0.666666666666667; mGzi6[1][1]= 0.166666666666667; mGzi6[1][2]= -0.577350269189626;
    mGzi6[2][0]= 0.166666666666667; mGzi6[2][1]= 0.666666666666667; mGzi6[2][2]= -0.577350269189626;
    mGzi6[3][0]= 0.166666666666667; mGzi6[3][1]= 0.166666666666667; mGzi6[3][2]=  0.577350269189626;
    mGzi6[4][0]= 0.666666666666667; mGzi6[4][1]= 0.166666666666667; mGzi6[4][2]=  0.577350269189626;
    mGzi6[5][0]= 0.166666666666667; mGzi6[5][1]= 0.666666666666667; mGzi6[5][2]=  0.577350269189626;
    // 9点 積分
    mGzi9[0][0]= 0.166666666666667; mGzi9[0][1]= 0.166666666666667; mGzi9[0][2]= -0.774596669241483;
    mGzi9[1][0]= 0.666666666666667; mGzi9[1][1]= 0.166666666666667; mGzi9[1][2]= -0.774596669241483;
    mGzi9[2][0]= 0.166666666666667; mGzi9[2][1]= 0.666666666666667; mGzi9[2][2]= -0.774596669241483;
    mGzi9[3][0]= 0.166666666666667; mGzi9[3][1]= 0.166666666666667; mGzi9[3][2]= 0.0;
    mGzi9[4][0]= 0.666666666666667; mGzi9[4][1]= 0.166666666666667; mGzi9[4][2]= 0.0;
    mGzi9[5][0]= 0.166666666666667; mGzi9[5][1]= 0.666666666666667; mGzi9[5][2]= 0.0;
    mGzi9[6][0]= 0.166666666666667; mGzi9[6][1]= 0.166666666666667; mGzi9[6][2]= 0.774596669241483;
    mGzi9[7][0]= 0.666666666666667; mGzi9[7][1]= 0.166666666666667; mGzi9[7][2]= 0.774596669241483;
    mGzi9[8][0]= 0.166666666666667; mGzi9[8][1]= 0.666666666666667; mGzi9[8][2]= 0.774596669241483;
    // 18点 積分
    mGzi18[0][0]=  0.166666666666667; mGzi18[0][1]=  0.166666666666667; mGzi18[0][2]= -0.774596669241483;
    mGzi18[1][0]=  0.166666666666667; mGzi18[1][1]=  0.666666666666667; mGzi18[1][2]= -0.774596669241483;          
    mGzi18[2][0]=  0.666666666666667; mGzi18[2][1]=  0.166666666666667; mGzi18[2][2]= -0.774596669241483;          
    mGzi18[3][0]=  0.000000000000000; mGzi18[3][1]=  0.500000000000000; mGzi18[3][2]= -0.774596669241483;          
    mGzi18[4][0]=  0.500000000000000; mGzi18[4][1]=  0.000000000000000; mGzi18[4][2]= -0.774596669241483;          
    mGzi18[5][0]=  0.500000000000000; mGzi18[5][1]=  0.500000000000000; mGzi18[5][2]= -0.774596669241483;          
    mGzi18[6][0]=  0.166666666666667; mGzi18[6][1]=  0.166666666666667; mGzi18[6][2]=  0.0;                         
    mGzi18[7][0]=  0.166666666666667; mGzi18[7][1]=  0.666666666666667; mGzi18[7][2]=  0.0;                         
    mGzi18[8][0]=  0.666666666666667; mGzi18[8][1]=  0.166666666666667; mGzi18[8][2]=  0.0;                         
    mGzi18[9][0]=  0.000000000000000; mGzi18[9][1]=  0.500000000000000; mGzi18[9][2]=  0.0;                         
    mGzi18[10][0]= 0.500000000000000; mGzi18[10][1]= 0.000000000000000; mGzi18[10][2]= 0.0;                         
    mGzi18[11][0]= 0.500000000000000; mGzi18[11][1]= 0.500000000000000; mGzi18[11][2]= 0.0;                         
    mGzi18[12][0]= 0.166666666666667; mGzi18[12][1]= 0.166666666666667; mGzi18[12][2]= 0.774596669241483;           
    mGzi18[13][0]= 0.166666666666667; mGzi18[13][1]= 0.666666666666667; mGzi18[13][2]= 0.774596669241483;           
    mGzi18[14][0]= 0.666666666666667; mGzi18[14][1]= 0.166666666666667; mGzi18[14][2]= 0.774596669241483;           
    mGzi18[15][0]= 0.000000000000000; mGzi18[15][1]= 0.500000000000000; mGzi18[15][2]= 0.774596669241483;           
    mGzi18[16][0]= 0.500000000000000; mGzi18[16][1]= 0.000000000000000; mGzi18[16][2]= 0.774596669241483;           
    mGzi18[17][0]= 0.500000000000000; mGzi18[17][1]= 0.500000000000000; mGzi18[17][2]= 0.774596669241483;

    //--
    //積分点の重み
    //--
    // 2点積分
    mW2[0]= 0.5; mW2[1]= 0.5;
    // 6点積分
    mW6[0]= 0.166666666666666; mW6[1]= 0.166666666666666; mW6[2]= 0.166666666666666;
    mW6[3]= 0.166666666666666; mW6[4]= 0.166666666666666; mW6[5]= 0.166666666666666;
    // 9点積分
    mW9[0]= 0.092592592592593; mW9[1]= 0.092592592592593; mW9[2]= 0.092592592592593;
    mW9[3]= 0.148148148148148; mW9[4]= 0.148148148148148; mW9[5]= 0.148148148148148;
    mW9[6]= 0.092592592592593; mW9[7]= 0.092592592592593; mW9[8]= 0.092592592592593;
    // 18点積分
    mW18[0]=  0.083333333333333; mW18[1]=  0.083333333333333; mW18[2]=  0.083333333333333;
    mW18[3]=  0.009259259259259; mW18[4]=  0.009259259259259; mW18[5]=  0.009259259259259;
    mW18[6]=  0.133333333333333; mW18[7]=  0.133333333333333; mW18[8]=  0.133333333333333;
    mW18[9]=  0.014814814814815; mW18[10]= 0.014814814814815; mW18[11]= 0.014814814814815;
    mW18[12]= 0.083333333333333; mW18[13]= 0.083333333333333; mW18[14]= 0.083333333333333;
    mW18[15]= 0.009259259259259; mW18[16]= 0.009259259259259; mW18[17]= 0.009259259259259;

    //積分点数
    mvIntegNum.resize(4);
    mvIntegNum[0]= 2;
    mvIntegNum[1]= 6;
    mvIntegNum[2]= 9;
    mvIntegNum[3]= 18;

    uiint numOfIntg,numOfShape,dof;
    //vector 領域確保
    //--
    //積分点2 節点6
    numOfIntg=2; numOfShape=6; dof=3;
    ResizeShape(mvN62, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr62, numOfIntg, numOfShape, dof);
    //積分点6 節点15
    numOfIntg=6; numOfShape=15; dof=3;
    ResizeShape(mvN156, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr156, numOfIntg, numOfShape, dof);
    //積分点9 節点15
    numOfIntg=9; numOfShape=15; dof=3;
    ResizeShape(mvN159, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr159, numOfIntg, numOfShape, dof);
    //積分点18 節点15
    numOfIntg=18; numOfShape=15; dof=3;
    ResizeShape(mvN1518, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr1518, numOfIntg, numOfShape, dof);

    //ヤコビアン
    mpJacobi62= new CJacobian;
    mpJacobi156= new CJacobian;
    mpJacobi159= new CJacobian;
    mpJacobi1518= new CJacobian;

    mpJacobi62->setupRegion(2,6);
    mpJacobi156->setupRegion(6,15);
    mpJacobi159->setupRegion(9,15);
    mpJacobi1518->setupRegion(18,15);
    //--
    //空間座標(x,y,z)の導関数配列確保
    //--
    // 6節点要素 積分点2
    numOfIntg= 2; numOfShape= 6; dof= 3;
    ResizeDeriv(mvdNdx62, numOfIntg, numOfShape, dof);
    // 15節点要素 積分点6
    numOfIntg= 6; numOfShape= 15; dof= 3;
    ResizeDeriv(mvdNdx156, numOfIntg, numOfShape, dof);
    // 15節点要素 積分点9
    numOfIntg= 9; numOfShape= 15; dof= 3;
    ResizeDeriv(mvdNdx159, numOfIntg, numOfShape, dof);
    // 15節点要素 積分点18
    numOfIntg= 18; numOfShape= 15; dof= 3;
    ResizeDeriv(mvdNdx1518, numOfIntg, numOfShape, dof);

    // detJの領域確保
    numOfIntg= 2; numOfShape= 6;
    Resize_detJ(mv_detJ62, numOfIntg);
    numOfIntg= 6; numOfShape= 15;
    Resize_detJ(mv_detJ156, numOfIntg);
    numOfIntg= 9; numOfShape= 15;
    Resize_detJ(mv_detJ159, numOfIntg);
    numOfIntg= 18; numOfShape= 15;
    Resize_detJ(mv_detJ1518, numOfIntg);
    
    //--
    //形状関数:N[積分点][i] : i番の形状関数
    //--
    setupShapeFunction(mvN62, 2, 6, mGzi2);     // 6節点   2点積分
    setupShapeFunction(mvN156, 6, 15, mGzi6);   // 15節点  6点積分
    setupShapeFunction(mvN159, 9, 15, mGzi9);   // 15節点  9点積分
    setupShapeFunction(mvN1518, 18, 15, mGzi18);// 15節点 18点積分
    

    //--
    //導関数  :dNdr[積分点][i][積分座標]
    //--
    setupShapeDeriv(mvdNdr62, 2, 6, mGzi2);// 6節点 2点積分
    setupShapeDeriv(mvdNdr156, 6, 15, mGzi6);// 15節点 6点積分
    setupShapeDeriv(mvdNdr159, 9, 15, mGzi9);// 15節点 9点積分
    setupShapeDeriv(mvdNdr1518, 18, 15, mGzi18);// 15節点 18点積分


    //形状関数積分値 セットアップ
    mvIntegValue15.resize(15);
    mvIntegValue6.resize(6);
    setupIntegValue15();
    setupIntegValue6();
}
CShapePrism::~CShapePrism()
{
    ;
}


// 形状関数のセットアップ
//
void CShapePrism::setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3])
{
    uiint igauss;
    double L1,L2,ze;

    switch(numOfShape){
        case(6):// 6節点
            for(igauss=0; igauss< numOfIntg; igauss++){
                L1= Gzi[igauss][0]; L2= Gzi[igauss][1]; ze= Gzi[igauss][2];

                ShapeFunction6(N, igauss, L1,L2,ze);
            };
            break;
        
        case(15):// 15節点
            for(igauss=0; igauss< numOfIntg; igauss++){
                L1= Gzi[igauss][0]; L2= Gzi[igauss][1]; ze= Gzi[igauss][2];

                ShapeFunction15(N, igauss, L1,L2,ze);
            };
            break;
            
        default:
            break;
    }
}

// 導関数のセットアップ
//
void CShapePrism::setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3])
{
    uiint igauss;
    double L1,L2,ze;

    switch(numOfShape){
        case(6):
            for(igauss=0; igauss< numOfIntg; igauss++){
                L1= Gzi[igauss][0]; L2= Gzi[igauss][1];  ze= Gzi[igauss][2];

                ShapeDeriv6(dNdr,igauss, L1,L2,ze);
            };
            break;

        case(15):
            for(igauss=0; igauss< numOfIntg; igauss++){
                L1= Gzi[igauss][0]; L2= Gzi[igauss][1];  ze= Gzi[igauss][2];

                ShapeDeriv15(dNdr,igauss, L1,L2,ze);
            };
            break;

        default:
            break;
    }
}


//形状関数
//--
void CShapePrism::ShapeFunction6(vvdouble& N, const uiint& igauss,
                 const double& L1, const double& L2, const double& gzi)
{
    double a= 1.0 - L1 - L2;
    N[igauss][0]= 0.50* a*(1.0 - gzi);
    N[igauss][1]= 0.50*L1*(1.0 - gzi);
    N[igauss][2]= 0.50*L2*(1.0 - gzi);
    N[igauss][3]= 0.50* a*(1.0 + gzi);
    N[igauss][4]= 0.50*L1*(1.0 + gzi);
    N[igauss][5]= 0.50*L2*(1.0 + gzi);
}
void CShapePrism::ShapeFunction15(vvdouble& N, const uiint& igauss,
                 const double& L1, const double& L2, const double& gzi)
{
    double a= 1.0 - L1 - L2;

    N[igauss][0]=  0.5* a*(1.0-gzi)*(2.0*a -2.0-gzi);
    N[igauss][1]=  0.5*L1*(1.0-gzi)*(2.0*L1-2.0-gzi);
    N[igauss][2]=  0.5*L2*(1.0-gzi)*(2.0*L2-2.0-gzi);
    N[igauss][3]=  0.5* a*(1.0+gzi)*(2.0*a -2.0+gzi);
    N[igauss][4]=  0.5*L1*(1.0+gzi)*(2.0*L1-2.0+gzi);
    N[igauss][5]=  0.5*L2*(1.0+gzi)*(2.0*L2-2.0+gzi);
    N[igauss][6]=  2.0*L1*a*(1.0-gzi);
    N[igauss][7]=  2.0*L1*L2*(1.0-gzi);
    N[igauss][8]=  2.0*L2*a*(1.0-gzi);
    N[igauss][9]=  2.0*L1*a*(1.0+gzi);
    N[igauss][10]= 2.0*L1*L2*(1.0+gzi);
    N[igauss][11]= 2.0*L2*a*(1.0+gzi);
    N[igauss][12]=  a*(1.0-gzi*gzi);
    N[igauss][13]= L1*(1.0-gzi*gzi);
    N[igauss][14]= L2*(1.0-gzi*gzi);
}


//導関数
//--
void CShapePrism::ShapeDeriv6(vvvdouble& dNdr, const uiint& igauss,
                 const double& L1, const double& L2, const double& gzi)
{
    double a= 1.0 - L1 - L2;

    // dN/dr
    //
    dNdr[igauss][0][0]= -0.50*(1.0-gzi);
    dNdr[igauss][1][0]=  0.50*(1.0-gzi);
    dNdr[igauss][2][0]=  0.0;
    dNdr[igauss][3][0]= -0.50*(1.0+gzi);
    dNdr[igauss][4][0]=  0.50*(1.0+gzi);
    dNdr[igauss][5][0]=  0.0;

    // dN/ds
    //
    dNdr[igauss][0][1]= -0.50*(1.0-gzi);
    dNdr[igauss][1][1]=  0.0;
    dNdr[igauss][2][1]=  0.50*(1.0-gzi);
    dNdr[igauss][3][1]= -0.50*(1.0+gzi);
    dNdr[igauss][4][1]=  0.0;
    dNdr[igauss][5][1]=  0.50*(1.0+gzi);

    // dN/dt
    //
    dNdr[igauss][0][2]= -0.50*a;
    dNdr[igauss][1][2]= -0.50*L1;
    dNdr[igauss][2][2]= -0.50*L2;
    dNdr[igauss][3][2]=  0.50*a;
    dNdr[igauss][4][2]=  0.50*L1;
    dNdr[igauss][5][2]=  0.50*L2;
}
void CShapePrism::ShapeDeriv15(vvvdouble& dNdr, const uiint& igauss,
                 const double& L1, const double& L2, const double& gzi)
{
    double a= 1.0- L1 - L2;

    // dN/dr
    //
    dNdr[igauss][0][0]=  -0.5*(1.0-gzi)*(4.0*a -gzi-2.0);
    dNdr[igauss][1][0]=   0.5*(1.0-gzi)*(4.0*L1-gzi-2.0);
    dNdr[igauss][2][0]=   0.0;
    dNdr[igauss][3][0]=  -0.5*(1.0+gzi)*(4.0*a +gzi-2.0);
    dNdr[igauss][4][0]=   0.5*(1.0+gzi)*(4.0*L1+gzi-2.0);
    dNdr[igauss][5][0]=   0.0;
    dNdr[igauss][6][0]=   2.0*(1.0-gzi)*(1.0-2.0*L1-L2);
    dNdr[igauss][7][0]=   2.0*L2*(1.0-gzi);
    dNdr[igauss][8][0]=  -2.0*L2*(1.0-gzi);
    dNdr[igauss][9][0]=   2.0*(1.0+gzi)*(1.0-2.0*L1-L2);
    dNdr[igauss][10][0]=  2.0*L2*(1.0+gzi);
    dNdr[igauss][11][0]=  -2.0*L2*(1.0+gzi);
    dNdr[igauss][12][0]=  -(1.0-gzi*gzi);
    dNdr[igauss][13][0]=   (1.0-gzi*gzi);
    dNdr[igauss][14][0]=   0.0;

    // dN/ds
    //
    dNdr[igauss][0][1]= -0.5*(1.0-gzi)*(4.0*a -gzi-2.0);
    dNdr[igauss][1][1]=  0.0;
    dNdr[igauss][2][1]=  0.5*(1.0-gzi)*(4.0*L2-gzi-2.0);
    dNdr[igauss][3][1]= -0.5*(1.0+gzi)*(4.0*a +gzi-2.0);
    dNdr[igauss][4][1]=  0.0;
    dNdr[igauss][5][1]=  0.5*(1.0+gzi)*(4.0*L2+gzi-2.0);
    dNdr[igauss][6][1]= -2.0*L1*(1.0-gzi);
    dNdr[igauss][7][1]=  2.0*L1*(1.0-gzi);
    dNdr[igauss][8][1]=  2.0*(1.0-gzi)*(1.0-L1-2.0*L2);
    dNdr[igauss][9][1]= -2.0*L1*(1.0+gzi);
    dNdr[igauss][10][1]=  2.0*L1*(1.0+gzi);
    dNdr[igauss][11][1]=  2.0*(1.0+gzi)*(1.0-L1-2.0*L2);
    dNdr[igauss][12][1]= -(1.0-gzi*gzi);
    dNdr[igauss][13][1]=  0.00;
    dNdr[igauss][14][1]=  (1.0-gzi*gzi);
    
    // dN/dt
    //
    dNdr[igauss][0][2]=   a*(L1+L2+gzi-0.5);
    dNdr[igauss][1][2]=  L1*(-L1+gzi+0.5);
    dNdr[igauss][2][2]=  L2*(-L2+gzi+0.5);
    dNdr[igauss][3][2]=   a*(-L1-L2+gzi+0.5);
    dNdr[igauss][4][2]=  L1*(L1+gzi-0.5);
    dNdr[igauss][5][2]=  L2*(L2+gzi-0.5);
    dNdr[igauss][6][2]= -2*L1*a;
    dNdr[igauss][7][2]= -2*L1*L2;
    dNdr[igauss][8][2]= -2*L2*a;
    dNdr[igauss][9][2]=   2*L1*a;
    dNdr[igauss][10][2]=  2*L1*L2;
    dNdr[igauss][11][2]=  2*L2*a;
    dNdr[igauss][12][2]= -2*a*gzi;
    dNdr[igauss][13][2]= -2*L1*gzi;
    dNdr[igauss][14][2]= -2*L2*gzi;
}


// 返り値:重みの配列, 引数:積分点数
// ----
double* CShapePrism::Weight(const uiint& integNum)
{
    switch(integNum){
        case(2):
            return mW2;
            break;
        case(6):
            return mW6;
            break;
        case(9):
            return mW9;
            break;
        case(18):
            return mW18;
            break;
        default:
            
            break;
    }
}
// 返り値:重み, 引数:重みの配列Index
// ----
double& CShapePrism::Weight_pt2(const uiint& igauss){ return mW2[igauss];}
double& CShapePrism::Weight_pt6(const uiint& igauss){ return mW6[igauss];}
double& CShapePrism::Weight_pt9(const uiint& igauss){  return mW9[igauss];}
double& CShapePrism::Weight_pt18(const uiint& igauss){ return mW18[igauss];}

// 形状関数   引数:積分点位置index, 形状関数番号(節点)
// ----
double& CShapePrism::N62(const uiint& igauss, const uiint& ishape){ return mvN62[igauss][ishape];}
double& CShapePrism::N156(const uiint& igauss, const uiint& ishape){ return mvN156[igauss][ishape];}
double& CShapePrism::N159(const uiint& igauss, const uiint& ishape){  return mvN159[igauss][ishape];}
double& CShapePrism::N1518(const uiint& igauss, const uiint& ishape){  return mvN1518[igauss][ishape];}

// 導関数   引数:積分点位置index, 形状関数番号(節点),微分方向(axis:0,1,2)
// ----
double& CShapePrism::dNdr62(const uiint& igauss, const uiint& ishape, const uiint& axis){ return mvdNdr62[igauss][ishape][axis];}
double& CShapePrism::dNdr156(const uiint& igauss, const uiint& ishape, const uiint& axis){ return mvdNdr156[igauss][ishape][axis];}
double& CShapePrism::dNdr159(const uiint& igauss, const uiint& ishape, const uiint& axis){  return mvdNdr159[igauss][ishape][axis];}
double& CShapePrism::dNdr1518(const uiint& igauss, const uiint& ishape, const uiint& axis){  return mvdNdr1518[igauss][ishape][axis];}


// 6節点:Prism dNdxのセットアップ
// --
void CShapePrism::Calc_dNdx6(const uiint& numOfInteg, CElement *pElement)
{
    vector<CNode*> vNode= pElement->getNode();
    uiint numOfShape(6);
    uiint dof(3);//3次元: J(3*3)行列

    switch(numOfInteg){
        case(2):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr62, mpJacobi62, vNode, mvdNdx62, mv_detJ62);
            break;
        default:
            break;
    }
}
// 15節点:Prism2  dNdxのセットアップ
// --
void CShapePrism::Calc_dNdx15(const uiint& numOfInteg, CElement *pElement)
{
    // Elementの関数:辺ノードの取得
    //  CNode* getEdgeInterNode(const uint& edgeIndex){ return mvEdgeInterNode[edgeIndex];}

    vector<CNode*> vNode= pElement->getNode();
    uiint numOfShape(15);

//    CNode *pNode;
//    uint ishape,iedge;
//    // 要素の辺ノードをFrontISTR順に並び替えて,vNodeに追加
//    //
//    for(ishape=6; ishape< numOfShape; ishape++){
//        iedge= mpISTR2Edge->PrismEdgeNum(ishape);
//        pNode= pElement->getEdgeInterNode(iedge);
//        vNode.push_back(pNode);
//    };

    uiint dof(3);//3次元
    switch(numOfInteg){
        case(6):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr156, mpJacobi156, vNode, mvdNdx156, mv_detJ156);
            break;
        case(9):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr159, mpJacobi159, vNode, mvdNdx159, mv_detJ159);
            break;
        case(18):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr1518, mpJacobi1518, vNode, mvdNdx1518, mv_detJ1518);
            break;
    }
}

// detJ
//
double& CShapePrism::detJ(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss)
{
    switch(elemType){
        case(ElementType::Prism):
            if(numOfInteg==2)
                return mv_detJ62[igauss];
            break;
        case(ElementType::Prism2):
            if(numOfInteg==6)
                return mv_detJ156[igauss];
            if(numOfInteg==9)
                return mv_detJ159[igauss];
            if(numOfInteg==18)
                return mv_detJ1518[igauss];
            break;
    }
}



// 形状関数の積分値セットアップ：コンストラクタからコール
// --
// 15節点 :9点積分を利用
// --
void CShapePrism::setupIntegValue15()
{
    uiint ishape;
    for(ishape=0; ishape < 15; ishape++) mvIntegValue15[ishape]= 0.0;

    uiint igauss;
    for(igauss=0; igauss < 9; igauss++){
        for(ishape=0; ishape < 15; ishape++){
            mvIntegValue15[ishape] += mvN159[igauss][ishape];
        };
    };

    //Normalize
    for(ishape=0; ishape < 15; ishape++) mvIntegValue15[ishape] /= 9.0;
}
// --
// 6節点
// --
void CShapePrism::setupIntegValue6()
{
    uiint ishape;
    for(ishape=0; ishape < 6; ishape++) mvIntegValue6[ishape]= 1.0/6.0; /* 0.16666666666666 */
}












