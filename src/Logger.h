/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Logger.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _LOGGER_H_3a27862d_3f3f_4845_8145_6216d7c2faeb
#define	_LOGGER_H_3a27862d_3f3f_4845_8145_6216d7c2faeb
#include <iostream>
#include <cstdio> 
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
    uint myRank;
    ofstream ofs;
    string msPathName;   
    string msLogFileBaseName;
    string msLogFileName;
    vuint mvOutputDevice;
    uint  mnCurrentState;
    string  msModeString;  
    string  msOutputString;
    char    mcFillStr;     
    string  msLoggerPref;  
    double  mdDummyValue;
    uint    muDummyValue;
    int     miDummyValue;
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
    void InfoDisplay();
    void setFileBaseName(const string& filebase);
    void initializeLogFile(const uint& rank);
    void finalizeLogFile();
    void setMode(const uint& mode);
    void setProperty(const uint& mode, const uint& outputDevice);
    void Monitor(const uint& mode, const uint& id, const double& value, const string& message);
    void Monitor(const uint& mode, const uint& id, const vdouble& vValue, const string& message);
    void Monitor(const uint& mode, const uint& id, const vint& vValue, const string& message);
    void Info(const uint& mode, const string& message);
    void Info(const uint& mode, const string& message1, const string& message2);
    void Info(const uint& mode, const string& message, const uint& num);
    void Info(const uint& mode, const string& message, const double& val);
    double& getDDummyValue();
    uint&   getUDummyValue();
    int&    getIDummyValue();
};
}
#endif	/* _LOGGER_H_ */
