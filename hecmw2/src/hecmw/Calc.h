/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Calc.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
 */

#ifndef CE3E2AF_026D_4cff_9038_A465D0ED626D
#define CE3E2AF_026D_4cff_9038_A465D0ED626D

//#include <cstdio>
#include <cstdlib>
#include <cctype>  //char型 判定
#ifdef MSVC
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <algorithm> // transform

#include <iostream>
#include <stack>
#include <string>
using namespace std;

//#include <boost/lexical_cast.hpp>
//using namespace boost;

#include "Poland.h"

namespace pmw
{
class CCalc
{
private:
    CCalc();
public:
    static CCalc* Instance() {
        static CCalc moCalc;
        return &moCalc;
    }
    virtual ~CCalc();

    // 計算前に要素中心の座標値 & 境界値をワーク変数にセット
    void setElementParam(const double& val, const double& x, const double& y, const double& z);

private:
    //--
    // ワーク変数
    //--
    bool   mbParam;   // 座標がセットされているか
    double mX, mY, mZ;// 要素中心の座標値のワーク変数
    double mBndVal;   // 境界値のワーク変数
    //--
    // 文字列の判定
    //--
    bool isNumber(string token);
    bool isCoord(string token, double& val);//-- (x,y,z)の場合:座標値に置き換え
    bool isBndValue(string token);          //-- (val)の判定
    bool isPI(string token);                //-- (pi)の判定
    bool isOperator(string token);          //-- (演算子)の判定

public:
    //逆ポーランド式の計算
    double Exec(CPoland *pNum);
};
}
#endif //include guard














