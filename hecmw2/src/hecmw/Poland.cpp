/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Poland.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
 */

#include "Poland.h"
using namespace pmw;


CPoland::CPoland()
{
    mPolandForm= "";
    mOrigForm= "";

    //--
    //演算子の優先順位:Low?High
    //--

    //単項演算子
    mmPrio["sin"]= 4;
    mmPrio["cos"]= 4;
    mmPrio["tan"]= 4;

    mmPrio["asin"]= 4;
    mmPrio["acos"]= 4;
    mmPrio["atan"]= 4;

    mmPrio["sinh"]= 4;
    mmPrio["cosh"]= 4;
    mmPrio["tanh"]= 4;

    mmPrio["log"]= 4;
    mmPrio["log10"]= 4;

    mmPrio["sqrt"]= 4;

    //二項演算子
    mmPrio["+"]= 0;
    mmPrio["-"]= 1;
    mmPrio["*"]= 2;
    mmPrio["/"]= 3;

    mmPrio["^"]= 7;//--- pow

    // ( )
    mmPrio["("]= 8;
    mmPrio[")"]= 8;

    ////// 元 : "[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
    ////m_strNumber="[0-9]*\\.\?[0-9]+([eE][-+]\?[0-9]+)\?";
}
CPoland::~CPoland()
{
    ;
}
// 0-9 || "."
bool CPoland::isNumber(char c)
{
    bool bCheck(false);

    if( isdigit(c) ) bCheck=true;
    if( c == '.' ) bCheck=true;

    return bCheck;
}
// "("
bool CPoland::isLBracket(char c)
{
    bool bCheck(false);

    if( c== '(' ) bCheck=true;

    return bCheck;
}
// ")"
bool CPoland::isRBracket(char c)
{
    bool bCheck(false);

    if( c== ')' ) bCheck=true;

    return bCheck;
}
// a-z
bool CPoland::isAlpha(char c)
{
    bool bCheck(false);

    if( isalpha(c) ) bCheck=true;

    return bCheck;
}
// ^,*,/,+,-
bool CPoland::isOperator(char c)
{
    bool bCheck(false);

    if( c=='^' ) bCheck=true;
    if( c=='*' ) bCheck=true;
    if( c=='/' ) bCheck=true;
    if( c=='+' ) bCheck=true;
    if( c=='-' ) bCheck=true;

    return bCheck;
}

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
//
// asin, acos, atan
string CPoland::getA(string& str, size_t& i)
{
    string swork;

    swork = str[i];
    i++;
    if(str[i]=='s') {     //asin
        swork += str[i];
        i++;
        if(str[i]=='i') {
            swork += str[i];
            i++;
        }
        if(str[i]=='n') {
            swork += str[i];
            i++;
        }
    } else if(str[i]=='c') { //acos
        swork += str[i];
        i++;
        if(str[i]=='o') {
            swork += str[i];
            i++;
        }
        if(str[i]=='s') {
            swork += str[i];
            i++;
        }
    } else if(str[i]=='t') { //atan
        swork += str[i];
        i++;
        if(str[i]=='a') {
            swork += str[i];
            i++;
        }
        if(str[i]=='n') {
            swork += str[i];
            i++;
        }
    }

    return swork;
}
// sin, sinh, sqrt
string CPoland::getS(string& str, size_t& i)
{
    string swork;

    swork = str[i];
    i++;
    if(str[i]=='i') { // sin, sinh
        swork += str[i];
        i++;
        if(str[i]=='n') {
            swork += str[i];
            i++;
        }
        if(str[i]=='h') {
            swork += str[i];
            i++;
        }
    } else if(str[i]=='q') { // sqrt
        swork += str[i];
        i++;
        if(str[i]=='r') {
            swork += str[i];
            i++;
        }
        if(str[i]=='t') {
            swork += str[i];
            i++;
        }
    }

    return swork;
}
// cos, cosh
string CPoland::getC(string& str, size_t& i)
{
    string swork;

    swork = str[i];
    i++;
    if(str[i]=='o') {
        swork += str[i];
        i++;
    }
    if(str[i]=='s') {
        swork += str[i];
        i++;
    }
    if(str[i]=='h') {
        swork += str[i];
        i++;
    }

    return swork;
}
// tan, tanh
string CPoland::getT(string& str, size_t& i)
{
    string swork;

    swork = str[i];
    i++;
    if(str[i]=='a') {
        swork += str[i];
        i++;
    }
    if(str[i]=='n') {
        swork += str[i];
        i++;
    }
    if(str[i]=='h') {
        swork += str[i];
        i++;
    }

    return swork;
}
// log, log10
string CPoland::getL(string& str, size_t& i)
{
    string swork;

    swork = str[i];
    i++;
    if(str[i]=='o') {
        swork += str[i];
        i++;
    }
    if(str[i]=='g') {
        swork += str[i];
        i++;
    }
    if(str[i]=='1') {
        swork += str[i];
        i++;
        if(str[i]=='0') {
            swork += str[i];
            i++;
        } else {
            cout << "Error ----- log10になってない." << endl;
            return "log";
        }
    }

    return swork;
}
// val
string CPoland::getV(string& str, size_t& i)
{
    string swork;

    swork = str[i];
    i++;
    if(str[i]=='a') {
        swork += str[i];
        i++;
    }
    if(str[i]=='l') {
        swork += str[i];
        i++;
    }

    return swork;
}
// pi
string CPoland::getP(string& str, size_t& i)
{
    string swork;

    swork = str[i];
    i++;
    if(str[i]=='i') {
        swork += str[i];
        i++;
    }
    return swork;
}
//
// ファイルから取得した式は、スペースが有ったり無かったりする ⇒ トークン分割
//
void CPoland::setOrigForm(string str)
{
    size_t nLength= str.length();

    // エラー処理
    if(nLength == 0) {
        cout << "Error ------ CPoland::setOrigForm" << endl;
        return;
    }

    // 全て小文字へ変換
    transform(str.begin(), str.end(), str.begin(), (int(*)(int))std::tolower);
    //transform(str.begin(), str.end(), str.begin(), ::tolower);

    mOrigForm.clear();

    //--
    // 先頭が単項演算子"-","+"の場合: 先頭に"0 "を挿入
    //--
    if(str[0]=='-' || str[0]=='+') mOrigForm= "0 ";

    ////    //--
    ////    // 正規表現:指数形式を含む数字を全て取得 :テスト
    ////    //--
    ////    size_t nCountSum = searchNumber(str);
    ////    size_t nCount=0;

    //--
    // 数字・文字列:カッコ( ),[ sqrt, sin, cos, tan, ... ], [ val, pi, x, y, z], [ ^, *, /, +, - ]
    //--
    size_t i=0;
    while(i < nLength) {

        string swork;
        // 数字
        if( isNumber(str[i]) ) {

            //指数形式を含む数字取得
            string sNum= getSciNumber(str,i);
            swork = sNum;

            ////    //テスト:正規表現結果との比較
            ////    if(nCountSum==nCount){
            ////        cout<< "error ---- Poland::setOrigForm, regex_result count do not match" <<endl;
            ////        break;
            ////    }
            ////    if(sNum==mv_regStr[nCount]){
            ////        cout << "Poland::setOrigForm, fix searchNumber:" << sNum << ", count:" << nCount << endl;
            ////    }else{
            ////        cout << "warn ---- Poland::setOrigForm, regex_result and func_result do not match" << endl;
            ////    }
            ////    nCount++;

            ////debug
            //cout << "Poland::setOrigForm, Number:" << swork << endl;
        }
        // 二項演算子:^, *, /, +, -
        else if( isOperator(str[i]) ) {
            swork = str[i];
            i++;
            ////debug
            //cout << "Poland::setOrigForm, BinaryOp:" << swork << endl;
        }
        // 括弧 ()
        else if( isLBracket(str[i]) || isRBracket(str[i]) ) {
            swork = str[i];
            i++;
            ////debug
            //cout << "Poland::setOrigForm, Bracket:" << swork << endl;
        }
        // 英字
        //----
        // a: asin, acos, atan
        // s: sin, sinh, sqrt
        // c: cos, cosh
        // t: tan, tanh
        // l: log, log10
        //----
        // v: val, p: pi
        //----
        // x,  y,  z,
        //----
        else if( isAlpha(str[i]) ) {
            switch(str[i]) {
            case('a'):
                swork = getA(str, i);
                break;
            case('s'):
                swork = getS(str, i);
                break;
            case('c'):
                swork = getC(str, i);
                break;
            case('t'):
                swork = getT(str, i);
                break;
            case('l'):
                swork = getL(str, i);
                break;
            case('v'):
                swork = getV(str, i);
                break;
            case('p'):
                swork = getP(str, i);
                break;
            case('x'):
                swork = str[i];
                i++;
                break;
            case('y'):
                swork = str[i];
                i++;
                break;
            case('z'):
                swork = str[i];
                i++;
                break;
            default:
                break;
            }//switch end

            ////debug
            //cout << "Poland::setOrigForm, UnaryOp || Param:" << swork << endl;
        }
        mOrigForm += swork + ' ';// 後ろに空白を付けて mOrigForm に付け足す.

    };//while(nLength) end

    ////debug
    //cout << "Poland::setOrigForm, mOrigForm: " << mOrigForm << endl;

    mOrigForm += '\0';// 終端ヌル
}

bool CPoland::isNumber(string& str)
{
    bool bCheck(false);

    bCheck = isNumber(str[0]);

    return bCheck;
}
bool CPoland::isLBracket(string& str)
{
    bool bCheck(false);

    bCheck = isLBracket(str[0]);

    return bCheck;
}
bool CPoland::isRBracket(string& str)
{
    bool bCheck(false);

    bCheck = isRBracket(str[0]);

    return bCheck;
}
//--
// 単項演算子: Unary operator
//--
bool CPoland::isUnaryOp(string& str)
{
    bool bCheck(false);

    if("sqrt"==str) bCheck=true;

    if("sin"==str) bCheck=true;
    if("cos"==str) bCheck=true;
    if("tan"==str) bCheck=true;

    if("asin"==str) bCheck=true;
    if("acos"==str) bCheck=true;
    if("atan"==str) bCheck=true;

    if("sinh"==str) bCheck=true;
    if("cosh"==str) bCheck=true;
    if("tanh"==str) bCheck=true;

    if("log"==str) bCheck=true;
    if("log10"==str) bCheck=true;


    return bCheck;
}
//--
// 二項演算子: Binary operator
//--
bool CPoland::isBinaryOp(string& str)
{
    bool bCheck(false);

    if("^"==str) bCheck=true;
    if("*"==str) bCheck=true;
    if("/"==str) bCheck=true;
    if("+"==str) bCheck=true;
    if("-"==str) bCheck=true;

    return bCheck;
}
//--
// 変数: val, x,y,z && pi
//--
bool CPoland::isParam(string& str)
{
    bool bCheck(false);

    if("val"==str) {
        bCheck=true;
        return bCheck;
    } else if("x"==str) {
        bCheck=true;
        return bCheck;
    } else if("y"==str) {
        bCheck=true;
        return bCheck;
    } else if("z"==str) {
        bCheck=true;
        return bCheck;
    } else if("pi"==str) {
        bCheck=true;
        return bCheck;
    } else {
        return bCheck;
    }
}

//--
// Scientific_Number(指数形式)数字を含む数字の取得
// # isNumber==true の場合に呼び出し
//--
string CPoland::getSciNumber(string str, size_t& i)
{
    string sNumber;
    size_t nEpos=0;       //'e'位置
    bool   bEflag(false); //
    size_t nPMpos=0;      //'+-'位置
    bool   bPMflag(false);//

    sNumber= str[i];//--- 先頭がisNumberの場合に呼ばれる関数なので数字
    i++;
    for(; i < str.length(); i++) {
        //Error処理
        if(bPMflag) {
            if( !isNumber(str[i]) && nPMpos+1==i) {
                cout << "Error ------------- Poland::getSciNumber " << endl;
                exit(0);//------------'+-'の後ろが数値でない. --- exit(0);
            }
        }
        //数値判定
        if( isNumber(str[i]) ) {
            sNumber += str[i];
        } else if( isE(str[i]) ) {
            sNumber += str[i];
            nEpos = i;//------- 'e'位置
            bEflag= true;
        } else if( isPlusMinus(str[i]) ) {
            if( bEflag && i==nEpos+1 ) {
                sNumber += str[i];
                nPMpos = i;
                bPMflag= true;
            } else {
                break;    //----------------------- break!
            }
        } else {
            break;//--------------------- break!
        }
    };

    return sNumber;
}
bool CPoland::isE(char c)
{
    if(c=='e') {
        return true;
    } else {
        return false;
    }
}
bool CPoland::isPlusMinus(char c)
{
    if(c=='+' || c=='-') {
        return true;
    } else {
        return false;
    }
}

//////--
////// 正規表現による数字取得
//////--
////size_t CPoland::searchNumber(string str)
////{
////    /*
////    sregex regNumber=sregex::compile(m_strNumber);
////    sregex_iterator iter( str.begin(), str.end(), regNumber);
////    sregex_iterator last;
////    while ( iter != last ) {
////        for ( size_t i = 0; i < iter->size(); ++i ) {
////            cout << " pos:" << iter->position(i);
////            cout << " len:" << iter->length(i);
////            cout << " str:" << iter->str(i);
////        }
////        cout << endl;
////        ++iter;//--------------------- iter インクリメント
////    }
////    */
////    sregex regNumber=sregex::compile(m_strNumber);
////    sregex_iterator iter( str.begin(), str.end(), regNumber);
////    sregex_iterator last;
////
////    for(; iter != last; iter++){
////        mv_regPos.push_back( iter->position() );
////        mv_regLen.push_back( iter->length() );
////        mv_regStr.push_back( iter->str() );
////
////        cout << " pos:" << iter->position();
////        cout << " len:" << iter->length();
////        cout << " str:" << iter->str();
////        cout << endl;
////    };
////
////    return mv_regStr.size();
////}




//---
// 逆ポーランド記法へ変換
//---
void CPoland::transForm()
{
    stack<string> stk; //スタック(LIFO):演算子を比較
    string      token; //トークン

    stringstream ss(mOrigForm);//--元の式:FIFO

    mPolandForm.clear();//逆ポーランド式クリア
    mPolandForm="";
    //--
    // 元の数式のトークンを一つずつ判定
    //--
    while(!ss.eof()) {

        ss >> token;

        //---- break
        if(token[0]=='\0') break;

        //--
        // トークン判定
        //--
        if( isNumber(token) ) {
            //数字はそのまま保存
            mPolandForm += token + " ";

            ////debug
            //cout << "Poland::transForm Number:" << token << endl;
        } else if( isParam(token) ) {
            //パラメータは数字の代用なので、そのまま保存
            mPolandForm += token + " ";

            ////debug
            //cout << "Poland::transForm Param:" << token << endl;
        } else if( isLBracket(token) ) {
            //左カッコ"(" ⇒ スタックに積む
            stk.push(token);

            ////debug
            //cout << "Poland::transForm LBracket:" << token << endl;
        } else if( isBinaryOp(token) ) {
            //二項演算子
            if( !stk.empty() ) {
                //--
                // 先に積まれた演算子がある場合は比較する
                //--
                string fOp= stk.top();//-- 先に積まれた演算子:first op
                if( mmPrio[fOp] >= mmPrio[token] ) {
                    ////debug
                    //cout << "先行演算子 高い&同じ,  先行演算子:" << fOp << " token:" << token << endl;
                    //bool bCheck(false); if(fOp==token) bCheck=true;

                    if( mmPrio[fOp] != mmPrio["("] ) {
                        mPolandForm += stk.top() + " ";//--- 先に積まれた演算子の方が優先度が高い(同じ)
                        stk.pop();//------------------------ 破棄
                        stk.push(token);//------------------ 優先度の低い演算子をスタックに積む
                    } else {
                        stk.push(token);//------------------ 先に"("が積まれている
                    }
                    ////debug : *==* の様子
                    //if(bCheck){
                    //	cout << " *==* : Form  " << mPolandForm << endl;
                    //	cout << " stk: ";
                    //	for(size_t i=0; i < stk.size(); i++) cout << " " << stk.top();
                    //	cout << endl;
                    //}
                } else {
                    stk.push(token);//------------------ 先に積まれた演算子の優先度低い:スタックに積む

                    ////debug
                    //cout << "先に積まれた演算子の優先度低 " << fOp << " token:" << token << endl;

                }//if(演算子比較)end

            } else {
                ////debug
                //cout << "最初の演算子: " << token << endl;
                //--
                // 最初の演算子 ⇒ スタックに積む
                //--
                stk.push(token);
            }
            ////debug
            //cout << "Poland::transForm BinaryOp:" << token << endl;
        } else if( isUnaryOp(token)  ) {
            //単項演算子
            if( !stk.empty() ) {
                //--
                // 先に積まれた演算子がある場合は比較する
                //--
                string fOp= stk.top();//-- 先に積まれた演算子:first op
                if( mmPrio[fOp] >= mmPrio[token] ) {

                    if( mmPrio[fOp] != mmPrio["("] ) {
                        mPolandForm += stk.top() + " ";//--- 先に積まれた演算子の方が優先度が高い(同じ)
                        stk.pop();//------------------------ 破棄
                        stk.push(token);//------------------ 優先度の低い演算子をスタックに積む
                    } else {
                        stk.push(token);//------------------ 先に"("が積まれている
                    }

                } else {
                    stk.push(token);
                }//if(演算子比較)end
            } else {
                //--
                // 最初の演算子 ⇒ スタックに積む
                //--
                stk.push(token);
            }
            ////debug
            //cout << "Poland::transForm UnaryOp:" << token << endl;
        } else if( isRBracket(token) ) {
            ////debug
            //cout << "Poland::transForm RBracket:" << token << endl;

            // 右カッコ")" ⇒ 左カッコ"("までスタックから演算子を取得
            // ## 比較は既に実行されている
            while( !stk.empty() ) {
                if( stk.top()=="(" ) {
                    stk.pop();//破棄
                    break;
                } else {
                    ////debug
                    //cout << "LBracket==( までスタックから取得 演算子: " << stk.top() << endl;

                    mPolandForm += stk.top() + " ";
                    stk.pop();//破棄
                }
            };//while(!stk.empty) end
        } else {
            // Error
            cout << "Warn ------ Poland::transForm : 特定できないトークン : " << token << endl;
        }
    };//while(!ss.eof) end

    //--
    // 残りの演算子をFormに付け加える
    //--
    //cout << "残りの演算子 数:" << stk.size() << endl;
    //
    while( !stk.empty() ) {
        ////debug
        //cout << "残り演算子: " << stk.top() << endl;

        mPolandForm += stk.top() + " ";
        stk.pop();//破棄
    };

    mPolandForm += '\0';//終端ヌル
}

//--
// Debug: PolandFormの文字列出力
//--
string CPoland::outputForm()
{
    return mPolandForm;
}

//-----------------
// Calc(計算)利用関数 : 数式(逆ポーランド記法)を提供
//-----------------
// 1.Calc利用準備(SS化)
//--
void CPoland::prepForm()
{
    mssForm.str("");
    mssForm.clear();

    mssForm.str(mPolandForm);//---- Calc利用準備(SS化): preparation
}
//--
// 2.数式の終端判定
//-
bool CPoland::isEnd()
{
    return mssForm.eof();
}
//--
// 3.数式からトークン取得
//--
string CPoland::getToken()
{
    string token;

    mssForm >> token;

    return token;
}

