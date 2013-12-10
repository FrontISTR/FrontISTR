/*
 ----------------------------------------------------------
|
| Software Name : conv2mw3 Ver 0.1 beta
|
|   main.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include <string>
#include <iostream>
#include "ConvMain.h"

using namespace std;

int main(int argc, char* argv[])
{
        CConvMain* pConv= CConvMain::Instance();
        CMessage* pMsg= CMessage::Instance();

        string s_iName, s_oName, sOpt;
        bool   bOpFlag(false);

        if(argc <= 4 && argc > 2){
                s_iName = argv[1];// 入力ファイル名
                s_oName = argv[2];// 出力ファイル名

                if(argc==4){//オプション
                        sOpt= argv[3];
                        if( string::npos != sOpt.find("-u") && sOpt.length()==2){
                                bOpFlag=true;
                        }else{
                                pMsg->usage("conv2mw3 <input_filename> <output_filename> <[-u]>");//-u:ユーザー入力あり
                                exit(0);//-------------------------exit(0)
                        }
                }
        }else{
                pMsg->usage("conv2mw3 <input_filename> <output_filename> <[-u]>");//-u:ユーザー入力あり
                exit(0);//-------------------------exit(0)
        }

        pMsg->banner();
        //--
        // メッシュ・ファイル変換
        //--
        pConv->FileRead_FISTR4(s_iName, bOpFlag);// 入力:FrontISTR v4 メッシュ
        pConv->setupAssyModel();                 // データセットアップ(出力準備)
        pConv->FileWrite_MW3(s_oName);           // 出力:MW3 メッシュ

        return 1;
}
