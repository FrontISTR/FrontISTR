/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Logger.h
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
#ifndef _LOGGER_H_3a27862d_3f3f_4845_8145_6216d7c2faeb
#define	_LOGGER_H_3a27862d_3f3f_4845_8145_6216d7c2faeb
#include <iostream>
#include "CommonFile.h"
#include "CommonStd.h"
#include "TypeDef.h"
#include "LoggerType.h"
#include "Info.h"
#include <cstdarg> // ...
#include <cstring> //strlen
using namespace std;

namespace Utility
{
class CLogger
{
private:
    CLogger();
    uiint myRank;
    ofstream ofs;
    string msPathName;
    string msLogFileBaseName;
    string msLogFileName;
    vuint mvOutputDevice;
    uiint  mnCurrentState;
    string  msModeString;
    string  msOutputString;
    char    mcFillStr;
    string  msLoggerPref;
    double  mdDummyValue;
    uiint   muDummyValue;
    iint    miDummyValue;
    bool    BooleanMode(const uiint& mode);
    string& ModeToString(const uiint& mode);
    string& ArgToString(const uiint& id, const double& value);
    string& ArgToString(const uiint& id, const vdouble& vValue);
    string& ArgToString(const uiint& id, const vint& vValue);
    string& InfoToString(const uiint& num);
    string& InfoToString(const double& val);
    short mnWidth;
public:
    static CLogger* Instance() {
        static CLogger logger;
        return &logger;
    }
    virtual ~CLogger();
public:
    void InfoDisplay();

    void setFileBaseName(const string& filebase);

    void initializeLogFile(const uiint& rank);
    void finalizeLogFile();

    void setMode(const uiint& mode);
    void setProperty(const uiint& mode, const uiint& outputDevice);

    void Monitor(const uiint& mode, const uiint& id, const double& value, const string& message);
    void Monitor(const uiint& mode, const uiint& id, const vdouble& vValue, const string& message);
    void Monitor(const uiint& mode, const uiint& id, const vint& vValue, const string& message);

    void Info(const uiint& mode, const string& message);
    void Info(const uiint& mode, const string& message1, const string& message2);
    void Info(const uiint& mode, const string& message, const uiint& num);
    void Info(const uiint& mode, vstring& vMessage, vuint& vNumber);
    void Info(const uiint& mode, const string& message, const double& val);
    void Info(const uiint& mode, string format, vector<void*>& param );
    void Info_format(const uiint& mode, const char* format, ... );


    void setDDummyValue(const double& val);
    void setUDummyValue(const uiint& val);
    void setIDummyValue(const iint& val);

    double& getDDummyValue();
    uiint&  getUDummyValue();
    iint&   getIDummyValue();
};
}
#endif	/* _LOGGER_H_ */
