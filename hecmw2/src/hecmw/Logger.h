/* 
 * File:   Logger.h
 * Author: ktakeda
 *
 * Modify     2009/05/07
 * Created on 2009/05/01, 16:31
 */

#ifndef _LOGGER_H_3a27862d_3f3f_4845_8145_6216d7c2faeb
#define	_LOGGER_H_3a27862d_3f3f_4845_8145_6216d7c2faeb

#include <iostream>
#include <cstdio> //sprintf
//#include "boost/lexical_cast.hpp"

#include "CommonFile.h"
#include "CommonStd.h"
#include "TypeDef.h"

#include "LoggerType.h"
#include "Info.h"

namespace Utility{
class CLogger{
private:
    CLogger();
    
    uiint myRank;//自分のrank

    ofstream ofs;

    string msPathName;   // Logファイルのパス名
    string msLogFileBaseName;//Logファイルのベース名
    string msLogFileName;// Logファイル名

    ////map<string, string, less<string> > msOutputType; //Mode別の出力先(Monitor, Disk)

    vuint mvOutputDevice;// OutputDevice on Mode
    uiint  mnCurrentState;// Mode


    string  msModeString;  // Mode 表示用(display) string
    string  msOutputString;// Message表示用(display) string
    char    mcFillStr;     // 埋め文字
    string  msLoggerPref;  // Logger Errorの表示用接頭辞

    double  mdDummyValue;// 参照変数ダミー :Loggerを引用している各メソッドで参照変数ダミーが必要になった場合に使用
    uiint   muDummyValue;// 参照変数ダミー :
    iint    miDummyValue;// 参照変数ダミー :

    bool    BooleanMode(const uiint& mode);
    string& ModeToString(const uiint& mode);

    string& ArgToString(const uiint& id, const double& value);
    string& ArgToString(const uiint& id, const vdouble& vValue);
    string& ArgToString(const uiint& id, const vint& vValue);

    string& InfoToString(const uiint& num);
    string& InfoToString(const double& val);

    short mnWidth;

public:
    static CLogger* Instance(){
	static CLogger logger;
	return &logger;
    }
    virtual ~CLogger();

public:
    // Infomation Display
    void InfoDisplay();
    
    //// ランクの設定 => initializeLogFile
    //void setRank(const uint& rank){ myRank= rank;}

    // LogFile
    void setFileBaseName(const string& filebase);//Logファイル名の変更
    void initializeLogFile(const uiint& rank);
    void finalizeLogFile();

    // Logger Prop
    void setMode(const uiint& mode);
    void setProperty(const uiint& mode, const uiint& outputDevice);

    // Logger Method
    void Monitor(const uiint& mode, const uiint& id, const double& value, const string& message);
    void Monitor(const uiint& mode, const uiint& id, const vdouble& vValue, const string& message);
    void Monitor(const uiint& mode, const uiint& id, const vint& vValue, const string& message);
    void Info(const uiint& mode, const string& message);
    void Info(const uiint& mode, const string& message1, const string& message2);
    void Info(const uiint& mode, const string& message, const uiint& num);
    void Info(const uiint& mode, vstring& vMessage, vuint& vNumber);
    void Info(const uiint& mode, const string& message, const double& val);

    // Loggerを引用している各メソッドで参照変数ダミーが必要になった場合に使用
    void setDDummyValue(const double& val);
    void setUDummyValue(const uiint& val);
    void setIDummyValue(const iint& val);
    double& getDDummyValue();
    uiint&  getUDummyValue();
    iint&   getIDummyValue();
};
}
#endif	/* _LOGGER_H_ */

