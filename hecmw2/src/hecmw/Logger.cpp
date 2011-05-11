//
// Logger.cpp
//                  2009.05.08
//                  2009.05.01
//                  k.Takeda
#include "Logger.h"
using namespace Utility;

CLogger::CLogger()
{
    // パス名
    // systemから*.exeパスを取得する必要があるが,まだしてない.
    msPathName = "./";

    // LogFile デフォルトのファイルベース名
    msLogFileBaseName = "hec_mw3";
    
 
    // Mode デフォルト
    mnCurrentState = LoggerMode::Info;

    // OutputDevice デフォルト
    mvOutputDevice.resize(5);
    mvOutputDevice[LoggerMode::MWDebug]= LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Debug]  = LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Error]  = LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Warn]   = LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Info]   = LoggerDevice::Display;

    mcFillStr = '_';
    msLoggerPref  = "Logger    ";
    

    // Loggerを引用している各メソッドで参照変数ダミーが必要になった場合に使用
    mdDummyValue= 0.0;
    muDummyValue= 0;
    miDummyValue= 0;

    mnWidth = 66;//文字列出力幅
}

CLogger::~CLogger()
{
}


// Infomation Display Header
//
void CLogger::InfoDisplay()
{
    cout << Info::Header() << endl;
}


// Logファイルベース名の変更
//
void CLogger::setFileBaseName(const string& filebase)
{
    msLogFileBaseName = filebase;
}

// LogFile Open
//
void CLogger::initializeLogFile(const uiint& rank)
{
    myRank= rank;

    // LogFile名を作成
    msLogFileName = msPathName + msLogFileBaseName + ".";
    stringstream ss;
    ss << myRank;
    msLogFileName += ss.str();
    msLogFileName += ".log";

    // 古いLogを除去
    msOutputString = "Remove old LogFile";
    if(remove(msLogFileName.c_str()) == 0){
        //モード：MWDebugの場合 リムーブしたことを表示
        if(mnCurrentState==LoggerMode::MWDebug)
            cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
    }
    
    
    ofs.open(msLogFileName.c_str(), ios::out | ios::app);
}

// LogFile Close
//
void CLogger::finalizeLogFile()
{
    ofs.close();
}

//
// カレントModeの設定
//
void CLogger::setMode(const uiint& mode)
{
    switch(mode){
        case(LoggerMode::MWDebug):
            mnCurrentState = mode;
            break;
        case(LoggerMode::Debug):
            mnCurrentState = mode;
            break;
        case(LoggerMode::Error):
            mnCurrentState = mode;
            break;
        case(LoggerMode::Warn):
            mnCurrentState = mode;
            break;
        case(LoggerMode::Info):
            mnCurrentState = mode;
            break;
        default:
            // 引数のチェック
            cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << "Logger Mode mismatch, at setMode()" << endl;
            break;
    }
}

// StateType別に画面端末か、ディスクのどちらかに出力する -> プロパティの設定
//
//
void CLogger::setProperty(const uiint& mode, const uiint& outputDevice)
{
    // arg Check.
    if(outputDevice > LoggerDevice::Display){
        cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << "Logger Device mismatch, at setProperty()" << endl;
        return;
    }

    // Device Property setup, each Mode
    switch(mode){
        case(LoggerMode::MWDebug):
            mvOutputDevice[LoggerMode::MWDebug] = outputDevice;
            break;
        case(LoggerMode::Debug):
            mvOutputDevice[LoggerMode::Debug] = outputDevice;
            break;
        case(LoggerMode::Error):
            mvOutputDevice[LoggerMode::Error] = outputDevice;
            break;
        case(LoggerMode::Warn):
            mvOutputDevice[LoggerMode::Warn] = outputDevice;
            break;
        case(LoggerMode::Info):
            mvOutputDevice[LoggerMode::Info] = outputDevice;
            break;
        default:
            // arg Check.
            cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << "Logger Mode mismatch, at setProperty()" << endl;
            break;
    }    
}

// Mode範囲チェック
//
bool CLogger::BooleanMode(const uiint& mode)
{
    bool bCheck(false);

    if(mode >= LoggerMode::Invalid){
        cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << "Logger Mode mismatch, ";
        bCheck=false;
    }else{
        bCheck=true;
    }

    return bCheck;
}

// Mode(unsigned int) => string
//
string& CLogger::ModeToString(const uiint& mode)
{
    switch(mode){
        case(LoggerMode::MWDebug):
            msModeString = "MW_Debug  ";
            break;
        case(LoggerMode::Debug):
            msModeString = "Debug     ";
            break;
        case(LoggerMode::Error):
            msModeString = "Error     ";
            break;
        case(LoggerMode::Warn):
            msModeString = "Warn      ";
            break;
        case(LoggerMode::Info):
            msModeString = "Info      ";
            break;
        default:
            msModeString = "LogMode Mismatch";
            break;
    }
    return msModeString;
}

// 指定Modeに合わせて、メッセージ,整数,実数 を出力
//
// Monitor 1
//
void CLogger::Monitor(const uiint& mode, const uiint& id, const double& value, const string& message)
{
    // Mode Check.
    //
    if(!BooleanMode(mode)){
        cout <<  "at Monitor 1" << endl;
        return;
    }

    //sprintf(mcOutputString, "%u %e", id, value);
    //msOutputString = mcOutputString;
    //msOutputString += "@" + message;
    msOutputString = " ID=";
    ArgToString(id,value);
    msOutputString += " " + message;
    
    // 指定Modeでの動作の可否
    //
    if(mode >= mnCurrentState){
        // Terminal Display
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        // Disk File
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
// Monitor 2
//
void CLogger::Monitor(const uiint& mode, const uiint& id, const vdouble& vValue, const string& message)
{
    // Mode Check.
    //
    if(!BooleanMode(mode)){
        cout <<  "at Monitor 2" << endl;
        return;
    }

    msOutputString = " ID=";
    ArgToString(id,vValue);
    msOutputString += " " + message;

    // 指定Modeでの動作の可否
    //
    if(mode >= mnCurrentState){
        // Terminal Display
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        // Disk File
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
// Monitor 3
//
void CLogger::Monitor(const uiint& mode, const uiint& id, const vint& vValue, const string& message)
{
    // Mode Check.
    //
    if(!BooleanMode(mode)){
        cout <<  "at Monitor 3" << endl;
        return;
    }

    msOutputString = " ID=";
    ArgToString(id,vValue);
    msOutputString += " " + message;
    
    // 指定Modeでの動作の可否
    //
    if(mode >= mnCurrentState){
        uiint i;
        // Terminal Display
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        // Disk File
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}

//  lexical_cast による引数の文字列への変換
//  -----
//  type.1
//
string& CLogger::ArgToString(const uiint& id, const double& value)
{
    stringstream ss;
    ss << id;
    msOutputString += ss.str();
    msOutputString += ",";
    ss.clear(); ss.str("");
    ss << value;
    msOutputString += ss.str() + ",";
    
    return msOutputString;
}
//  type.2
//
string& CLogger::ArgToString(const uiint& id, const vdouble& vValue)
{
    stringstream ss;
    ss << id;
    msOutputString += ss.str();
    msOutputString += ",";
    
    for(uiint i=0; i < vValue.size(); i++){
        ss.clear(); ss.str("");
        ss << vValue[i];
        msOutputString += ss.str() + ",";
    };
    return msOutputString;
}
//  type.3
//
string& CLogger::ArgToString(const uiint& id, const vint& vValue)
{
    stringstream ss;
    ss << id;
    msOutputString += ss.str();
    msOutputString += ",";
    for(uiint i=0; i < vValue.size(); i++){
        ss.clear(); ss.str("");
        ss << vValue[i];
        msOutputString += ss.str() + ",";
    };
    return msOutputString;
}


// Info()用 文字列変換
//
string& CLogger::InfoToString(const uiint& num)
{
    stringstream ss;
    ss << num;
    msOutputString += ss.str();
    return msOutputString;
}
string& CLogger::InfoToString(const double& val)
{
    stringstream ss;
    ss << val;
    msOutputString += ss.str();
    return msOutputString;
}

// 指定Modeにあわせて、メッセージを出力
//
void CLogger::Info(const uiint& mode, const string& message)
{
    // Mode Check.
    //
    if(!BooleanMode(mode)){
        cout <<  "at Info()" << endl;
        return;
    }

    msOutputString = message;

    // 指定Modeでの動作可否判断
    //
    if(mode >= mnCurrentState){
        // Terminal Display
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        // Disk File
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            //debug
            //cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << "Disk write at CLogger::Info()" << endl;

            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}

// Info string型
//
void CLogger::Info(const uiint& mode, const string& message1, const string& message2)
{
    // Mode Check.
    //
    if(!BooleanMode(mode)){
        cout <<  "at Info()" << endl;
        return;
    }

    msOutputString = message1 + message2;

    // 指定Modeでの動作可否判断
    //
    if(mode >= mnCurrentState){
        // Terminal Display
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        // Disk File
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            //debug
            //cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << "Disk write at CLogger::Info()" << endl;

            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}

// Info unsigned int 型
//
void CLogger::Info(const uiint& mode, const string& message, const uiint& num)
{
    // Mode Check.
    //
    if(!BooleanMode(mode)){
        cout <<  "at Info()" << endl;
        return;
    }

    msOutputString = message;// + ", ";
    InfoToString(num);

    // 指定Modeでの動作可否判断
    //
    if(mode >= mnCurrentState){
        // Terminal Display
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        // Disk File
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            //debug
            //cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << "Disk write at CLogger::Info()" << endl;

            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
// Info vstring vuint 型
//
void CLogger::Info(const uiint& mode, vstring& vMessage, vuint& vNumber)
{
    // Mode Check.
    //
    if(!BooleanMode(mode)){
        cout <<  "at Info()" << endl;
        return;
    }

    // msOutputStringに蓄積
    msOutputString.clear();
    for(uiint i=0; i < vMessage.size(); i++){
        msOutputString += vMessage[i];
        InfoToString(vNumber[i]);
    };

    // 指定Modeでの動作可否判断
    //
    if(mode >= mnCurrentState){
        // Terminal Display
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        // Disk File
        if(LoggerDevice::Disk == mvOutputDevice[mode]){

            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}

// Info double 型
//
void CLogger::Info(const uiint& mode, const string& message, const double& val)
{
    // Mode Check.
    //
    if(!BooleanMode(mode)){
        cout <<  "at Info()" << endl;
        return;
    }

    msOutputString = message;// + ", ";
    InfoToString(val);

    // 指定Modeでの動作可否判断
    //
    if(mode >= mnCurrentState){
        // Terminal Display
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        // Disk File
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            //debug
            //cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << "Disk write at CLogger::Info()" << endl;

            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}


// Loggerを引用している各メソッドで参照変数ダミーが必要になった場合に使用
//
double& CLogger::getDDummyValue()
{
    return mdDummyValue;
}
uiint& CLogger::getUDummyValue()
{
    return muDummyValue;
}
iint& CLogger::getIDummyValue()
{
    return miDummyValue;
}
void CLogger::setDDummyValue(const double& val)
{
    mdDummyValue= val;
}
void CLogger::setUDummyValue(const uiint& val)
{
    muDummyValue= val;
}
void CLogger::setIDummyValue(const iint& val)
{
    miDummyValue= val;
}


