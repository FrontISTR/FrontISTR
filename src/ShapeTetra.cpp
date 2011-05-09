//
// ShapeTetra.cpp
//
//
//
//
//              2010.01.28
//              k.Takeda
#include "ShapeTetra.h"
using namespace pmw;

CShapeTetra::CShapeTetra()
{
//    // MW3 辺番号 -> FrontISTR 2次節点番号
//    mvEdge_2_ISTR[0]= 6;
//    mvEdge_2_ISTR[1]= 4;
//    mvEdge_2_ISTR[2]= 5;
//    mvEdge_2_ISTR[3]= 7;
//    mvEdge_2_ISTR[4]= 8;
//    mvEdge_2_ISTR[5]= 9;
//
//    // FrontISTR -> MW3 辺番号
//    uint i;
//    for(i=0; i< 4; i++){
//        mvISTR_2_Edge[i]= i;
//    };
//    mvISTR_2_Edge[4]= 1;
//    mvISTR_2_Edge[5]= 2;
//    mvISTR_2_Edge[6]= 0;
//    mvISTR_2_Edge[7]= 3;
//    mvISTR_2_Edge[8]= 4;
//    mvISTR_2_Edge[9]= 5;
    
    //--
    //積分座標(体積座標)
    //--
    // 1点 積分
    // L1            :    L2            :   L3
    //-------------------------------------------
    mGzi1[0][0]= 0.25; mGzi1[0][1]= 0.25; mGzi1[0][2]= 0.25;
    
    // 4点 積分
    // L1                         :  L2                           :  L3
    //---------------------------------------------------------------------
    mGzi4[0][0]= 0.138196601125011; mGzi4[0][1]= 0.138196601125011; mGzi4[0][2]= 0.138196601125011;
    mGzi4[1][0]= 0.585410196624968; mGzi4[1][1]= 0.138196601125011; mGzi4[1][2]= 0.138196601125011;
    mGzi4[2][0]= 0.138196601125011; mGzi4[2][1]= 0.585410196624968; mGzi4[2][2]= 0.138196601125011;
    mGzi4[3][0]= 0.138196601125011; mGzi4[3][1]= 0.138196601125011; mGzi4[3][2]= 0.585410196624968;
    
    //15点 積分
    // L1                           : L2                              : L3
    //-----------------------------------------------------------------------
    mGzi15[0][0]=  0.25;              mGzi15[0][1]=  0.25;              mGzi15[0][2]=  0.25;
    mGzi15[1][0]=  0.091971078052723; mGzi15[1][1]=  0.091971078052723; mGzi15[1][2]=  0.091971078052723;
    mGzi15[2][0]=  0.724086765841831; mGzi15[2][1]=  0.091971078052723; mGzi15[2][2]=  0.091971078052723;
    mGzi15[3][0]=  0.091971078052723; mGzi15[3][1]=  0.724086765841831; mGzi15[3][2]=  0.091971078052723;
    mGzi15[4][0]=  0.091971078052723; mGzi15[4][1]=  0.091971078052723; mGzi15[4][2]=  0.724086765841831;
    mGzi15[5][0]=  0.319793627829630; mGzi15[5][1]=  0.319793627829630; mGzi15[5][2]=  0.319793627829630;
    mGzi15[6][0]=  0.040619116511110; mGzi15[6][1]=  0.319793627829630; mGzi15[6][2]=  0.319793627829630;
    mGzi15[7][0]=  0.319793627829630; mGzi15[7][1]=  0.040619116511110; mGzi15[7][2]=  0.319793627829630;
    mGzi15[8][0]=  0.319793627829630; mGzi15[8][1]=  0.319793627829630; mGzi15[8][2]=  0.040619116511110;
    mGzi15[9][0]=  0.056350832689629; mGzi15[9][1]=  0.056350832689629; mGzi15[9][2]=  0.443649167310371;
    mGzi15[10][0]= 0.443649167310371; mGzi15[10][1]= 0.056350832689629; mGzi15[10][2]= 0.056350832689629;
    mGzi15[11][0]= 0.443649167310371; mGzi15[11][1]= 0.443649167310371; mGzi15[11][2]= 0.056350832689629;
    mGzi15[12][0]= 0.056350832689629; mGzi15[12][1]= 0.443649167310371; mGzi15[12][2]= 0.443649167310371;
    mGzi15[13][0]= 0.056350832689629; mGzi15[13][1]= 0.443649167310371; mGzi15[13][2]= 0.056350832689629;
    mGzi15[14][0]= 0.443649167310371; mGzi15[14][1]= 0.056350832689629; mGzi15[14][2]= 0.443649167310371;

    //--
    // 重み
    //--
    // 1点 積分
    mW1[0]= 0.166666666666667;
    // 4点 積分
    mW4[0]= 0.041666666666667;
    mW4[1]= 0.041666666666667;
    mW4[2]= 0.041666666666667;
    mW4[3]= 0.041666666666667;
    // 15点 積分
    mW15[0]=  0.019753086419753; mW15[1]=  0.011989513963170; mW15[2]=  0.011989513963170;
    mW15[3]=  0.011989513963170; mW15[4]=  0.011989513963170; mW15[5]=  0.011511367871045;
    mW15[6]=  0.011511367871045; mW15[7]=  0.011511367871045; mW15[8]=  0.011511367871045;
    mW15[9]=  0.008818342151675; mW15[10]= 0.008818342151675; mW15[11]= 0.008818342151675;
    mW15[12]= 0.008818342151675; mW15[13]= 0.008818342151675; mW15[14]= 0.008818342151675;

    //--
    //積分点数
    //--
    mvIntegNum.resize(3);
    mvIntegNum[0]= 1;
    mvIntegNum[1]= 4;
    mvIntegNum[2]= 15;

    uint numOfIntg,numOfShape,dof;
    //--
    //vector配列確保
    //--
    // 4節点要素  積分点1
    numOfIntg=1; numOfShape=4; dof=3;
    ResizeShape(mvN41, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr41, numOfIntg,numOfShape,dof);
    // 10節点要素  積分点1
    numOfIntg=1; numOfShape=10; dof=3;
    ResizeShape(mvN101, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr101, numOfIntg,numOfShape,dof);
    // 10節点要素  積分点4
    numOfIntg=4; numOfShape=10; dof=3;
    ResizeShape(mvN104, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr104, numOfIntg,numOfShape,dof);
    // 10節点要素 積分点15
    numOfIntg=15; numOfShape=10; dof=3;
    ResizeShape(mvN1015, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr1015, numOfIntg,numOfShape,dof);

    //ヤコビアン
    mpJacobi41= new CJacobian;
    mpJacobi101= new CJacobian;
    mpJacobi104= new CJacobian;
    mpJacobi1015= new CJacobian;

    mpJacobi41->setupRegion(1,4);
    mpJacobi101->setupRegion(1,10);
    mpJacobi104->setupRegion(4,10);
    mpJacobi1015->setupRegion(15,10);
    //--
    //空間座標(x,y,z)の導関数配列確保
    //--
    // 4節点要素 積分点1
    numOfIntg= 1; numOfShape= 4; dof= 3;
    ResizeDeriv(mvdNdx41, numOfIntg, numOfShape, dof);
    // 10節点要素 積分点1
    numOfIntg= 1; numOfShape= 10; dof= 3;
    ResizeDeriv(mvdNdx101, numOfIntg, numOfShape, dof);
    // 10節点要素 積分点4
    numOfIntg= 4; numOfShape= 10; dof= 3;
    ResizeDeriv(mvdNdx104, numOfIntg, numOfShape, dof);
    // 10節点要素 積分点15
    numOfIntg= 15; numOfShape= 10; dof= 3;
    ResizeDeriv(mvdNdx1015, numOfIntg, numOfShape, dof);

    // detJの領域確保
    numOfIntg= 1; numOfShape= 4;
    Resize_detJ(mv_detJ41, numOfIntg);
    numOfIntg= 1; numOfShape= 10;
    Resize_detJ(mv_detJ101, numOfIntg);
    numOfIntg= 4; numOfShape= 10;
    Resize_detJ(mv_detJ104, numOfIntg);
    numOfIntg= 15; numOfShape= 10;
    Resize_detJ(mv_detJ1015, numOfIntg);


    //--
    //形状関数:N[積分点][i] : i番の形状関数
    //--
    setupShapeFunction(mvN41, 1,4,mGzi1);// 4節点 1点積分
    setupShapeFunction(mvN101, 1,10,mGzi1);// 10節点 1点積分
    setupShapeFunction(mvN104, 4,10,mGzi4);// 10節点 4点積分
    setupShapeFunction(mvN1015, 15,10,mGzi15);// 10節点 15点積分
    
    //--
    //導関数  :dNdr[積分点][i][vol_coord]
    //--
    uint igauss= 0;
    {
        ShapeDeriv4(mvdNdr41,igauss);// 4節点 1点積分
    }
    setupShapeDeriv(mvdNdr101, 1,10,mGzi1);// 10節点 1点積分
    setupShapeDeriv(mvdNdr104, 4,10,mGzi4);// 10節点 4点積分
    setupShapeDeriv(mvdNdr1015, 15,10,mGzi15);// 10節点 15点積分



    //--
    //形状関数積分値のセットアップ
    //--
    mvIntegValue10.resize(10);
    mvIntegValue4.resize(4);
    setupIntegValue10();
    setupIntegValue4();
}
CShapeTetra::~CShapeTetra()
{
    ;
}


//形状関数 & 導関数のセットアップ
//
void CShapeTetra::setupShapeFunction(vvdouble& N, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][3])
{
    double L1,L2,L3;
    uint igauss;
    
    switch(numOfShape){
        case(4):
            for(igauss=0; igauss< numOfIntg; igauss++){
                L1= Gzi[igauss][0]; L2= Gzi[igauss][1]; L3= Gzi[igauss][2];

                ShapeFunction4(N,igauss, L1,L2,L3);
            };
            break;
        
        case(10):
            for(igauss=0; igauss< numOfIntg; igauss++){
                L1= Gzi[igauss][0]; L2= Gzi[igauss][1]; L3= Gzi[igauss][2];

                ShapeFunction10(N,igauss, L1,L2,L3);
            };
            break;
        
        default:
            break;
    }
}
void CShapeTetra::setupShapeDeriv(vvvdouble& dNdr, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][3])
{
    double L1,L2,L3;
    uint igauss;
    
    switch(numOfShape){
        case(4)://4節点
            for(igauss=0; igauss< numOfIntg; igauss++){
                L1= Gzi[igauss][0]; L2= Gzi[igauss][1];  L3= Gzi[igauss][2];

                ShapeDeriv4(dNdr,igauss);
            };
            break;
            
        case(10)://10節点
            for(igauss=0; igauss< numOfIntg; igauss++){
                L1= Gzi[igauss][0]; L2= Gzi[igauss][1];  L3= Gzi[igauss][2];

                ShapeDeriv10(dNdr,igauss, L1,L2,L3);
            };
            break;
            
        default:
            break;
    }
}


// 形状関数 4 Node
//
void CShapeTetra::ShapeFunction4(vvdouble& N, const uint& igauss,
        const double& L1, const double& L2, const double& L3)
{
    N[igauss][0]= 1.0 - L1 - L2 - L3;
    N[igauss][1]= L1;
    N[igauss][2]= L2;
    N[igauss][3]= L3;
}
// 形状関数 10 Node
//
void CShapeTetra::ShapeFunction10(vvdouble& N, const uint& igauss, 
        const double& L1, const double& L2, const double& L3)
{
    double a= 1.0 - L1 - L2 - L3;

    N[igauss][0]= (2.0*a - 1.0)*a;
    N[igauss][1]= L1*(2.0*L1 - 1.0);
    N[igauss][2]= L2*(2.0*L2 - 1.0);
    N[igauss][3]= L3*(2.0*L3 - 1.0);
    N[igauss][4]= 4.0*L1*a;
    N[igauss][5]= 4.0*L1*L2;
    N[igauss][6]= 4.0*L2*a;
    N[igauss][7]= 4.0*L3*a;
    N[igauss][8]= 4.0*L1*L3;
    N[igauss][9]= 4.0*L2*L3;
}


// 導関数 4 Node
//
void CShapeTetra::ShapeDeriv4(vvvdouble& dNdr, const uint& igauss)
{
    dNdr[igauss][0][0] = -1.0;
    dNdr[igauss][1][0] =  1.0;
    dNdr[igauss][2][0] =  0.0;
    dNdr[igauss][3][0] =  0.0;

    dNdr[igauss][0][1] = -1.0;
    dNdr[igauss][1][1] =  0.0;
    dNdr[igauss][2][1] =  1.0;
    dNdr[igauss][3][1] =  0.0;

    dNdr[igauss][0][2] = -1.0;
    dNdr[igauss][1][2] =  0.0;
    dNdr[igauss][2][2] =  0.0;
    dNdr[igauss][3][2] =  1.0;
}
// 導関数 10 Node
//
void CShapeTetra::ShapeDeriv10(vvvdouble& dNdr, const uint& igauss,
        const double& L1, const double& L2, const double& L3)
{
    double a= 1.0 - L1 - L2 - L3;

    // dN/dL1
    dNdr[igauss][0][0]=  1.0 - 4.0*a;
    dNdr[igauss][1][0]=  4.0*L1 - 1.0;
    dNdr[igauss][2][0]=  0.0;
    dNdr[igauss][3][0]=  0.0;
    dNdr[igauss][4][0]=  4.0*(1.0 - 2.0*L1 - L2 - L3);
    dNdr[igauss][5][0]=  4.0*L2;
    dNdr[igauss][6][0]= -4.0*L2;
    dNdr[igauss][7][0]= -4.0*L3;
    dNdr[igauss][8][0]=  4.0*L3;
    dNdr[igauss][9][0]=  0.0;

    // dN/dL2
    dNdr[igauss][0][1]=  1.0 - 4.0*a;
    dNdr[igauss][1][1]=  0.0;
    dNdr[igauss][2][1]=  4.0*L2-1.0;
    dNdr[igauss][3][1]=  0.0;
    dNdr[igauss][4][1]= -4.0*L1;
    dNdr[igauss][5][1]=  4.0*L1;
    dNdr[igauss][6][1]=  4.0*(1.0 - L1 - 2.0*L2 - L3);
    dNdr[igauss][7][1]= -4.0*L3;
    dNdr[igauss][8][1]=  0.0;
    dNdr[igauss][9][1]=  4.0*L3;

    // dN/dL3
    dNdr[igauss][0][2]=  1.0 - 4.0*a;
    dNdr[igauss][1][2]=  0.0;
    dNdr[igauss][2][2]=  0.0;
    dNdr[igauss][3][2]=  4.0*L3 - 1.0;
    dNdr[igauss][4][2]= -4.0*L1;
    dNdr[igauss][5][2]=  0.0;
    dNdr[igauss][6][2]= -4.0*L2;
    dNdr[igauss][7][2]=  4.0*(1.0- L1 - L2 - 2.0*L3);
    dNdr[igauss][8][2]=  4.0*L1;
    dNdr[igauss][9][2]=  4.0*L2;
}

// 形状関数   引数:積分点位置index, 形状関数番号(節点)
// --
double& CShapeTetra::N41(const uint& igauss, const uint& ishape){ return mvN41[igauss][ishape];}
double& CShapeTetra::N101(const uint& igauss, const uint& ishape){ return mvN101[igauss][ishape];}
double& CShapeTetra::N104(const uint& igauss, const uint& ishape){ return mvN104[igauss][ishape];}
double& CShapeTetra::N1015(const uint& igauss, const uint& ishape){ return mvN1015[igauss][ishape];}

// 導関数   引数:積分点位置index, 形状関数番号(節点),微分方向(axis:0,1,2)
// --
double& CShapeTetra::dNdr41(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdr41[igauss][ishape][axis];}
double& CShapeTetra::dNdr101(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdr101[igauss][ishape][axis];}
double& CShapeTetra::dNdr104(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdr104[igauss][ishape][axis];}
double& CShapeTetra::dNdr1015(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdr1015[igauss][ishape][axis];}

// 返り値:重みの配列, 引数:積分点数
// --
double* CShapeTetra::Weight(const uint& integNum)
{
    switch(integNum){
        case(1):
            return mW1;
            break;
        case(4):
            return mW4;
            break;
        case(15):
            return mW15;
            break;
        default:
            break;
    }
}
// 返り値:重み, 引数:重みの配列Index
double& CShapeTetra::Weight_pt1() { return mW1[0];}
double& CShapeTetra::Weight_pt4(const uint& igauss){ return mW4[igauss];}
double& CShapeTetra::Weight_pt15(const uint& igauss){ return mW15[igauss];}


// 4節点:Tetra dNdxのセットアップ
// --
void CShapeTetra::Calc_dNdx4(const uint& numOfInteg, CElement *pElement)
{
    vector<CNode*> vNode= pElement->getNode();
    uint numOfShape(4);
    uint dof(3);//3次元: J(3*3)行列

    switch(numOfInteg){
        case(1):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr41, mpJacobi41, vNode, mvdNdx41, mv_detJ41);
            break;
        default:
            break;
    }
}
// 10節点:Tetra2  dNdxのセットアップ
// --
void CShapeTetra::Calc_dNdx10(const uint& numOfInteg, CElement *pElement)
{
    // Elementの関数:辺ノードの取得
    //  CNode* getEdgeInterNode(const uint& edgeIndex){ return mvEdgeInterNode[edgeIndex];}

    vector<CNode*> vNode= pElement->getNode();

    CNode *pNode;
    uint ishape,iedge;
    uint numOfShape(10);
    // 要素の辺ノードをFrontISTR順に並び替えて,vNodeに追加
    //
    for(ishape=4; ishape< numOfShape; ishape++){
        iedge= mpISTR2Edge->TetraEdgeNum(ishape);
        pNode= pElement->getEdgeInterNode(iedge);
        vNode.push_back(pNode);
    };
    uint dof(3);//3次元
    switch(numOfInteg){
        case(1):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr101, mpJacobi101, vNode, mvdNdx101, mv_detJ101);
            break;
        case(4):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr104, mpJacobi104, vNode, mvdNdx104, mv_detJ104);
            break;
        case(15):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr1015, mpJacobi1015, vNode, mvdNdx1015, mv_detJ1015);
            break;
    }
}

// detJ
//
double& CShapeTetra::detJ(const uint& elemType, const uint& numOfInteg, const uint& igauss)
{
    switch(elemType){
        case(ElementType::Tetra):
            if(numOfInteg==1)
                return mv_detJ41[igauss];
            break;
        case(ElementType::Tetra2):
            if(numOfInteg==1)
                return mv_detJ101[igauss];
            if(numOfInteg==4)
                return mv_detJ104[igauss];
            if(numOfInteg==15)
                return mv_detJ1015[igauss];
            break;
    }
}



// 形状関数の積分値セットアップ：コンストラクタからコール
// --
// 10節点 :4点積分を利用
// --
void CShapeTetra::setupIntegValue10()
{
    uint ishape;
    for(ishape=0; ishape < 10; ishape++) mvIntegValue10[ishape]= 0.0;
    
    uint igauss;
    for(igauss=0; igauss < 4; igauss++){
        for(ishape=0; ishape < 10; ishape++){
            mvIntegValue10[ishape] += mvN104[igauss][ishape];
        };
    };

    //Normalize
    for(ishape=0; ishape < 4; ishape++) mvIntegValue10[ishape] /= 4.0;
}
// --
// 4節点
// --
void CShapeTetra::setupIntegValue4()
{
    uint ishape;
    for(ishape=0; ishape < 4; ishape++) mvIntegValue4[ishape]= 0.25;
}












