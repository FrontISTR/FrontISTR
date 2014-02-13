/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryNode.cpp
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
#include <vector>
#include "BoundaryNode.h"
using namespace pmw;
CBoundaryNode::CBoundaryNode()
{
    ;
}
CBoundaryNode::~CBoundaryNode()
{
    ;
}
////void CBoundaryNode::resizeValue(const uiint& nNumOfDiffLevel, vuint& vDOF)
////{
////    //--
////    //引数(nNumOfDiffLevel)はコースグリッドBNodeの場合 : (最大Level + 1) == Factory::mMGLevel + 1
////    //--
////    mvValue.resize(nNumOfDiffLevel);
////    mvEntValue.resize(nNumOfDiffLevel);//ディレクレ基礎データ(Ent:エンティティ)
////
////    mvMarkingValue.resize(nNumOfDiffLevel);//Rcap_Refineのディレクレ分配時に利用
////    //--
////    // 境界値の初期化(全Level域、全DOF)
////    //--
////    for(uiint i=0; i < nNumOfDiffLevel; i++){
////        for(uiint idof=0; idof < vDOF.size(); idof++){
////            mvValue[i][idof]=0.0;
////            mvEntValue[i][idof]=0.0;
////        }
////    }
////}
void CBoundaryNode::resizeValue(const uiint& nNumOfDiffLevel)
{
    //--
    //引数(nNumOfDiffLevel)はコースグリッドBNodeの場合 : (最大Level + 1) == Factory::mMGLevel + 1
    //--
    mvValue.resize(nNumOfDiffLevel);
    mvEntValue.resize(nNumOfDiffLevel);//ディレクレ基礎データ(Ent:エンティティ)

    mvMarkingValue.resize(nNumOfDiffLevel);//Rcap_Refineのディレクレ分配時に利用
}
//--
// Neumann値(ノイマン)のaddのための初期化(0.0)
//--
void CBoundaryNode::initValue(const uiint& dof, const uiint& mgLevel)
{
    uiint diffLevel= mgLevel - mMGLevel;//自身のLevelとの差

    mvValue[diffLevel][dof]= 0.0;
    mvEntValue[diffLevel][dof]= 0.0;//ディレクレ基礎データ(Ent:エンティティ)
}

//--
// Rcapのディレクレ処理のために利用:自身のLevel〜MaxLevelまでBOOL初期化
//--
void CBoundaryNode::initRcapBool(const uiint& dof, const uiint& mgMaxLevel)
{
    uiint diffMaxLevel= mgMaxLevel - mMGLevel;//自身Levelとの差

    for(uiint i=0; i < diffMaxLevel+1; i++) {
        mvMarkingValue[i][dof]= false;//--- Rcap_Refineのディレクレ分配時に利用
    };
}

//--
// 境界値
//--
void CBoundaryNode::setValue(const uiint& dof, const uiint& mgLevel, const double& val)
{
    uiint diffLevel= mgLevel - mMGLevel;

    mvValue[diffLevel][dof]= val;

    mvMarkingValue[diffLevel][dof]= true;//Rcap_Refineのディレクレ分配時に利用
}
void CBoundaryNode::addValue(const uiint& dof, const uiint& mgLevel, const double& val)
{
    uiint diffLevel= mgLevel - mMGLevel;

    mvValue[diffLevel][dof] += val;
}
double& CBoundaryNode::getValue(const uiint& dof, const uiint& mgLevel)
{
    uiint diffLevel= mgLevel - mMGLevel;

    return mvValue[diffLevel][dof];
}
//
// Rcap_Refineのディレクレ分配時に利用
//
bool CBoundaryNode::isSetupValue(const uiint& dof, const uiint& mgLevel)
{
    uiint diffLevel= mgLevel - mMGLevel;

    return mvMarkingValue[diffLevel][dof];
}

//--
// ディレクレ基礎データ(Ent:エンティティ)
//--
void CBoundaryNode::setEntValue(const uiint& dof, const uiint& mgLevel, const double& val)
{
    uiint diffLevel= mgLevel - mMGLevel;

    mvEntValue[diffLevel][dof]= val;//Dirichletの数式計算値(Ent:エンティティ)
}
double& CBoundaryNode::getEntValue(const uiint& dof, const uiint& mgLevel)
{
    uiint diffLevel= mgLevel - mMGLevel;

    return mvEntValue[diffLevel][dof];
}


double& CBoundaryNode::getX()
{
    return mpNode->getX();
}
double& CBoundaryNode::getY()
{
    return mpNode->getY();
}
double& CBoundaryNode::getZ()
{
    return mpNode->getZ();
}


