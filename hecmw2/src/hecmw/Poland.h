/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Poland.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
 */

#ifndef CD89AB8_B86E_4b94_91BE_DA0ADEF32142
#define CD89AB8_B86E_4b94_91BE_DA0ADEF32142

#include <iostream>
#include <sstream>
#include <string>
#include <cctype> //char判定

#include <map>
#include <vector>
#include <stack>
#include <algorithm>
using namespace std;

////#include <boost/xpressive/xpressive.hpp>
////using namespace boost::xpressive;

namespace pmw
{
class CPoland
{
public:
    CPoland();
    virtual ~CPoland();

private:
    string mOrigForm;  // 元の数式
    string mPolandForm;// 逆ポーランド記法 数式

    map<string, short> mmPrio;// 演算子の優先順位

    stringstream mssForm;

    //--
    // ファイルに記述された数式 ⇒ 元の数式(トークン分割された数式)
    //--
    bool isNumber(char c);  // 0-9 || "."
    bool isLBracket(char c);// (
    bool isRBracket(char c);// )
    bool isAlpha(char c);   // a-z
    bool isOperator(char c);// ^,*,/,+,-

    // ----
    // a: asin, acos, atan
    // s: sin, sinh, sqrt
    // c: cos, cosh
    // t: tan, tanh
    // l: log, log10
    //
    // v: val
    // p: pi
    // x: x
    // y: y
    // z: z
    // ----
    string getA(string& str, size_t& i);
    string getS(string& str, size_t& i);
    string getC(string& str, size_t& i);
    string getT(string& str, size_t& i);
    string getL(string& str, size_t& i);
    string getV(string& str, size_t& i);
    string getP(string& str, size_t& i);

    //--
    // 元の数式 ⇒ 逆ポーランド式
    //--
    bool isNumber(string& str);
    bool isLBracket(string& str);
    bool isRBracket(string& str);

    bool isUnaryOp(string& str); //--単項演算子: Unary operator
    bool isBinaryOp(string& str);//--二項演算子: Binary operator

    bool isParam(string& str);//-----変数:val,x,y,z && pi

    //--
    // Scientific_Number(指数形式)数字を含む数字の取得
    //--
    string getSciNumber(string str, size_t& i);
    bool   isE(char c);
    bool   isPlusMinus(char c);

    ////    //--
    ////    // 正規表現による数値の取得(テスト用)
    ////    //--
    ////    size_t searchNumber(string str);
    ////
    ////    string m_strNumber;      //正規表現文字列
    ////    vector<size_t> mv_regPos;//結果:位置
    ////    vector<size_t> mv_regLen;//結果:文字列長
    ////    vector<string> mv_regStr;//結果:数字


public:
    //---
    // 逆ポーランド記法へ変換
    //---
    void setOrigForm(string str);//元の数式をセット: Numerical Formula
    void transForm();            //逆ポーランド記法へ変換

    string outputForm();//---- PolandFormの文字列出力

    //---
    // Calc(計算)で利用 : 数式(逆ポーランド記法)を提供
    //---
    void prepForm();//---- Calc利用準備(SS化): preparation

    bool isEnd();//------- 数式の終端判定
    string getToken();//-- 数式からトークン取得
};
}
#endif //include guard







