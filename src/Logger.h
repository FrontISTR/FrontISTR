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
#include "boost/lexical_cast.hpp"

#include "CommonFile.h"
#include "CommonStd.h"
#include "TypeDef.h"

#include "LoggerType.h"
#include "Info.h"

namespace Utility{
class CLogger{
private:
    CLogger();
    
    uint myRank;//自分のrank

    ofstream ofs;

    string msPathName;   // Logファイルのパス名
    string msLogFileBaseName;//Logファイルのベース名
    string msLogFileName;// Logファイル名

    ////map<string, string, less<string> > msOutputType; //Mode別の出力先(Monitor, Disk)

    vuint mvOutputDevice;// OutputDevice on Mode
    uint  mnCurrentState;// Mode


    string  msModeString;  // Mode 表示用(display) string
    string  msOutputString;// Message表示用(display) string
    char    mcFillStr;// 埋め文字
    string  msLoggerPref; // Logger Errorの表示用接頭辞

    bool    BooleanMode(const uint& mode);
    string& ModeToString(const uint& mode);

    string& ArgToString(const uint& id, const double& value);
    string& ArgToString(const uint& id, const vdouble& vValue);
    string& ArgToString(const uint& id, const vint& vValue);

    string& InfoToString(const uint& num);
    string& InfoToString(const double& val);

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
    void initializeLogFile(const uint& rank);
    void finalizeLogFile();

    // Logger Prop
    void setMode(const uint& mode);
    void setProperty(const uint& mode, const uint& outputDevice);

    // Logger Method
    void Monitor(const uint& mode, const uint& id, const double& value, const string& message);
    void Monitor(const uint& mode, const uint& id, const vdouble& vValue, const string& message);
    void Monitor(const uint& mode, const uint& id, const vint& vValue, const string& message);
    void Info(const uint& mode, const string& message);
    void Info(const uint& mode, const string& message1, const string& message2);
    void Info(const uint& mode, const string& message, const uint& num);
    void Info(const uint& mode, const string& message, const double& val);
};
}
#endif	/* _LOGGER_H_ */

