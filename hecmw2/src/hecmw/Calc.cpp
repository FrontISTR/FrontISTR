/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Calc.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
 */

#include "Calc.h"
using namespace pmw;

CCalc::CCalc()
{
    mX=0.0;
    mY=0.0;
    mZ=0.0;

    mBndVal=0.0;

    mbParam=false;
}
CCalc::~CCalc()
{
    ;
}
bool CCalc::isNumber(string token)
{
    bool bCheck(false);

    if( isdigit(token[0]) ) {
        bCheck= true;
        return bCheck;
    }
    if( token[0]== '.' ) {
        bCheck= true;
        return bCheck;
    }

    return bCheck;
}
void CCalc::setElementParam(const double& val, const double& x, const double& y, const double& z)
{
    mBndVal= val;

    mX= x;
    mY= y;
    mZ= z;

    mbParam= true;//要素座標 & 境界値がセットされた
}
bool CCalc::isCoord(string token, double& val)
{
    bool bCheck(false);

    if( token[0]== 'x' || token[0]== 'X' ) {
        val= mX;
        bCheck= true;
        return bCheck;
    }
    if( token[0]== 'y' || token[0]== 'Y' ) {
        val= mY;
        bCheck= true;
        return bCheck;
    }
    if( token[0]== 'z' || token[0]== 'Z' ) {
        val= mZ;
        bCheck= true;
        return bCheck;
    }

    return bCheck;
}
bool CCalc::isBndValue(string token)
{
    bool bCheck(false);

    if( token=="val" || token=="VAL" ) {
        bCheck= true;
        return bCheck;
    }

    return bCheck;
}
bool CCalc::isPI(string token)
{
    bool bCheck(false);

    if( token=="PI" || token=="pi" ) {
        bCheck=true;
        return bCheck;
    } else {
        return bCheck;
    }
}
bool CCalc::isOperator(string token)
{
    bool bCheck(false);

    if("sqrt"==token) bCheck=true;
    if("sin"==token)  bCheck=true;
    if("cos"==token)  bCheck=true;
    if("tan"==token)  bCheck=true;
    if("asin"==token) bCheck=true;
    if("acos"==token) bCheck=true;
    if("atan"==token) bCheck=true;
    if("sinh"==token) bCheck=true;
    if("cosh"==token) bCheck=true;
    if("tanh"==token) bCheck=true;
    if("log"==token)  bCheck=true;
    if("log10"==token) bCheck=true;
    if("^"==token) bCheck=true;
    if("*"==token) bCheck=true;
    if("/"==token) bCheck=true;
    if("+"==token) bCheck=true;
    if("-"==token) bCheck=true;

    return bCheck;
}
//--
// 逆ポーランド式の計算 : !注)powの乗数は整数に変換する
//--
double CCalc::Exec(CPoland *pNum)
{
    //////debug : 呼び出し回数を見るため
    ////cout << "-- Calc::Exec --" << endl;

    // Error:要素座標がセットされていない
    if(!mbParam) {
        cout << "Error ------------------- Calc::Proc, not Element Coord setting" << endl;
        return 0.0;
    }

    string token;      // トークン:逆ポーランド式から取得するトークン
    stack<double> stk; // 演算のスタック
    double d1,d2;      // ポップする値
    int    n1;         // "d1"を整数に変換するときの変数(pow用途)

    pNum->prepForm();  // 逆ポーランド式のストリーム準備

    //--
    // 逆ポーランド式の終わりまでループ
    //--
    while( !pNum->isEnd() ) {

        token= pNum->getToken();// 逆ポーランド式からトークン取り出し.

        if(token[0]=='\0') break;// ss終端

        ////debug
        //cout << "Calc::Exec token:" << token << endl;


        double val_c;
        //--
        // トークン評価: 逆ポーランド式の文字列
        //--
        if( isNumber(token) ) {
            // 数字:スタックに値を保存(push)
            double val= atof( token.c_str() );
            stk.push( val );
        } else if( isCoord(token, val_c) ) {
            // 座標値"x,y,z"のどれか: スタックに値を保存(push)
            stk.push( val_c );
        } else if( isBndValue(token) ) {
            // 境界値: スタックに値を保存(push)
            stk.push( mBndVal );
        } else if( isPI(token) ) {
            // π: スタックに値を保存(push)
            stk.push( M_PI );
        } else if( isOperator(token) ) {
            // 演算子:スタックから値を取得(pop)して計算 ⇒ 計算結果をスタックに積む(push)
            if(stk.size() < 1 ) {
                cout << "Error   Calc::Exec  数式エラー" << endl;
                exit(-1);
            }

            d1 = stk.top();// 引数1
            stk.pop();

            // 単項演算
            if(token=="sqrt") {
                stk.push( sqrt(d1) );
            } else if(token=="sin") {
                stk.push( sin(d1) );
            } else if(token=="cos") {
                stk.push( cos(d1) );
            } else if(token=="tan") {
                stk.push( tan(d1) );
            } else if(token=="asin") {
                stk.push( asin(d1) );
            } else if(token=="acos") {
                stk.push( acos(d1) );
            } else if(token=="atan") {
                stk.push( atan(d1) );
            } else if(token=="sinh") {
                stk.push( sinh(d1) );
            } else if(token=="cosh") {
                stk.push( cosh(d1) );
            } else if(token=="tanh") {
                stk.push( tanh(d1) );
            } else if(token=="log" ) {
                stk.push( log(d1)  );
            } else if(token=="log10") {
                stk.push( log10(d1));
            }
            // 2項演算
            else {
                d2 = stk.top();// 引数2
                stk.pop();

                if(token=="+") {
                    stk.push( d2 + d1 );
                } else if(token=="-") {
                    stk.push( d2 - d1 );
                } else if(token=="*") {
                    stk.push( d2 * d1 );
                } else if(token=="/") {
                    stk.push( d2 / d1 );
                } else if(token=="^") {
                    n1=(int)d1;    //乗数は整数(int)に変換
                    stk.push( pow(d2, n1) );
                }
                // 定義されていない記号
                else {
                    cout << "Error  Calc::Exec  定義されていない記号: " << token << endl;
                    exit(-1);
                }
            }
        } else {
            // 不明 token
            cout << "Warn  Calc::Exec  不明 token:" << token << endl;

        }//if(数値)
    };//while end

    mbParam= false;//計算が終わったのでElementパラメータ値を無効にしておく.

    return stk.top();//---- 計算結果
}









