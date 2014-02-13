/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Logger.cpp
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
#include "Logger.h"
using namespace Utility;
CLogger::CLogger()
{
    msPathName = "./";
    msLogFileBaseName = "hec_mw3";

    mnCurrentState = LoggerMode::Info;
    mvOutputDevice.resize(5);
    mvOutputDevice[LoggerMode::MWDebug]= LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Debug]  = LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Error]  = LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Warn]   = LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Info]   = LoggerDevice::Display;

    mcFillStr = '_';
    msLoggerPref  = "Logger    ";

    mdDummyValue=  0.0;
    muDummyValue=  100;
    miDummyValue= -100;

    mnWidth = 66;
}
CLogger::~CLogger()
{
}
void CLogger::InfoDisplay()
{
    cout << Info::Header() << endl;
}
void CLogger::setFileBaseName(const string& filebase)
{
    msLogFileBaseName = filebase;
}
void CLogger::initializeLogFile(const uiint& rank)
{
    myRank= rank;
    msLogFileName = msPathName + msLogFileBaseName + ".";
    stringstream ss;
    ss << myRank;
    msLogFileName += ss.str();
    msLogFileName += ".log";
    msOutputString = "Remove old LogFile";
    if(remove(msLogFileName.c_str()) == 0) {
        if(mnCurrentState==LoggerMode::MWDebug)
            cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
    }
    ofs.open(msLogFileName.c_str(), ios::out | ios::app);
}
void CLogger::finalizeLogFile()
{
    ofs.close();
}
void CLogger::setMode(const uiint& mode)
{
    switch(mode) {
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
        cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << "Logger Mode mismatch, at setMode()" << endl;
        break;
    }
}
void CLogger::setProperty(const uiint& mode, const uiint& outputDevice)
{
    if(outputDevice > LoggerDevice::Display) {
        cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << "Logger Device mismatch, at setProperty()" << endl;
        return;
    }
    switch(mode) {
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
        cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << "Logger Mode mismatch, at setProperty()" << endl;
        break;
    }
}
bool CLogger::BooleanMode(const uiint& mode)
{
    bool bCheck(false);
    if(mode >= LoggerMode::Invalid) {
        cout << msLoggerPref << setfill(mcFillStr) << setw(mnWidth) << "Logger Mode mismatch, ";
        bCheck=false;
    } else {
        bCheck=true;
    }
    return bCheck;
}
string& CLogger::ModeToString(const uiint& mode)
{
    switch(mode) {
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
void CLogger::Monitor(const uiint& mode, const uiint& id, const double& value, const string& message)
{
    if(!BooleanMode(mode)) {
        cout <<  "at Monitor 1" << endl;
        return;
    }
    msOutputString = " ID=";
    ArgToString(id,value);
    msOutputString += " " + message;
    if(mode >= mnCurrentState) {
        if(LoggerDevice::Display == mvOutputDevice[mode]) {
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]) {
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
void CLogger::Monitor(const uiint& mode, const uiint& id, const vdouble& vValue, const string& message)
{
    if(!BooleanMode(mode)) {
        cout <<  "at Monitor 2" << endl;
        return;
    }
    msOutputString = " ID=";
    ArgToString(id,vValue);
    msOutputString += " " + message;
    if(mode >= mnCurrentState) {
        if(LoggerDevice::Display == mvOutputDevice[mode]) {
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]) {
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
void CLogger::Monitor(const uiint& mode, const uiint& id, const vint& vValue, const string& message)
{
    if(!BooleanMode(mode)) {
        cout <<  "at Monitor 3" << endl;
        return;
    }
    msOutputString = " ID=";
    ArgToString(id,vValue);
    msOutputString += " " + message;
    if(mode >= mnCurrentState) {
        uiint i;
        if(LoggerDevice::Display == mvOutputDevice[mode]) {
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]) {
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
string& CLogger::ArgToString(const uiint& id, const double& value)
{
    stringstream ss;
    ss << id;
    msOutputString += ss.str();
    msOutputString += ",";
    ss.clear();
    ss.str("");
    ss << value;
    msOutputString += ss.str() + ",";
    return msOutputString;
}
string& CLogger::ArgToString(const uiint& id, const vdouble& vValue)
{
    stringstream ss;
    ss << id;
    msOutputString += ss.str();
    msOutputString += ",";
    for(uiint i=0; i < vValue.size(); i++) {
        ss.clear();
        ss.str("");
        ss << vValue[i];
        msOutputString += ss.str() + ",";
    };
    return msOutputString;
}
string& CLogger::ArgToString(const uiint& id, const vint& vValue)
{
    stringstream ss;
    ss << id;
    msOutputString += ss.str();
    msOutputString += ",";
    for(uiint i=0; i < vValue.size(); i++) {
        ss.clear();
        ss.str("");
        ss << vValue[i];
        msOutputString += ss.str() + ",";
    };
    return msOutputString;
}
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
void CLogger::Info(const uiint& mode, const string& message)
{
    if(!BooleanMode(mode)) {
        cout <<  "at Info()" << endl;
        return;
    }
    msOutputString = message;
    if(mode >= mnCurrentState) {
        if(LoggerDevice::Display == mvOutputDevice[mode]) {
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]) {
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
void CLogger::Info(const uiint& mode, const string& message1, const string& message2)
{
    if(!BooleanMode(mode)) {
        cout <<  "at Info()" << endl;
        return;
    }
    msOutputString = message1 + message2;
    if(mode >= mnCurrentState) {
        if(LoggerDevice::Display == mvOutputDevice[mode]) {
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]) {
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
void CLogger::Info(const uiint& mode, const string& message, const uiint& num)
{
    if(!BooleanMode(mode)) {
        cout <<  "at Info()" << endl;
        return;
    }
    msOutputString = message;
    InfoToString(num);
    if(mode >= mnCurrentState) {
        if(LoggerDevice::Display == mvOutputDevice[mode]) {
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]) {
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
void CLogger::Info(const uiint& mode, vstring& vMessage, vuint& vNumber)
{
    if(!BooleanMode(mode)) {
        cout <<  "at Info()" << endl;
        return;
    }
    msOutputString.clear();
    for(uiint i=0; i < vMessage.size(); i++) {
        msOutputString += vMessage[i];
        InfoToString(vNumber[i]);
    };
    if(mode >= mnCurrentState) {
        if(LoggerDevice::Display == mvOutputDevice[mode]) {
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]) {
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
void CLogger::Info(const uiint& mode, const string& message, const double& val)
{
    if(!BooleanMode(mode)) {
        cout <<  "at Info()" << endl;
        return;
    }
    msOutputString = message;
    InfoToString(val);
    if(mode >= mnCurrentState) {
        if(LoggerDevice::Display == mvOutputDevice[mode]) {
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]) {
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(mnWidth) << msOutputString << endl;
        }
    }
}
void CLogger::Info(const uiint& mode, string format, vector<void*>& param )
{
    if(!BooleanMode(mode)) {
        cout <<  "at Info()" << endl;
        return;
    }
    // Mode表示
    if(mode >= mnCurrentState) {
        if(LoggerDevice::Display == mvOutputDevice[mode])
            cout << ModeToString(mode) << right;
        if(LoggerDevice::Disk == mvOutputDevice[mode])
            ofs << ModeToString(mode)  << right;
    }

    // Value
    if(mode >= mnCurrentState) {

        iint *nVal;
        uiint *uVal;
        double* dVal;
        string sVal;
        uiint nLen= format.length();

        uiint ipara=0;
        for(uiint i=0; i < nLen; i++) {
            if(format[i] == '%') {
                ++i;
                switch( format[i] ) {
                case('d'):
                    nVal= (iint*)param[ipara];
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << dec << *nVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << dec << *nVal << " ";
                    ipara++;
                    break;
                case('u'):
                    uVal= (uiint*)param[ipara];
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << dec << *uVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << dec << *uVal << " ";
                    ipara++;
                    break;
                case('f'):
                    dVal= (double*)param[ipara];
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << fixed << *dVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << fixed << *dVal << " ";
                    ipara++;
                    break;
                case('e'):
                    dVal= (double*)param[ipara];
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << scientific << *dVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << scientific << *dVal << " ";
                    ipara++;
                    break;
                case('s'):
                    sVal= (char*)param[ipara];
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << sVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << sVal << " ";
                    ipara++;
                    break;
                default:
                    break;
                }
            }
        };
        if(LoggerDevice::Display == mvOutputDevice[mode])
            cout << endl;
        if(LoggerDevice::Disk == mvOutputDevice[mode])
            ofs << endl;

    }//if ( nCurrentState ) end

}
void CLogger::Info_format(const uiint& mode, const char* format, ... )
{
    if(!BooleanMode(mode)) {
        cout <<  "at Info()" << endl;
        return;
    }
    // Mode表示
    if(mode >= mnCurrentState) {
        if(LoggerDevice::Display == mvOutputDevice[mode])
            cout << ModeToString(mode) << right;
        if(LoggerDevice::Disk == mvOutputDevice[mode])
            ofs << ModeToString(mode)  << right;
    }

    // Value 出力
    //
    if(mode >= mnCurrentState) {
        uiint nLength = strlen(format);
        va_list list;
        va_start( list, format);
        for(uiint i=0; i < nLength; i++) {
            if(format[i] == '%') {
                ++i;
                switch( format[i] ) {
                case('d'): {
                    iint nVal = va_arg(list, iint);
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << dec << nVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << dec << nVal << " ";
                }
                break;
                case('u'): {
                    uiint nVal = va_arg(list, uiint);
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << dec << nVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << dec << nVal << " ";
                }
                break;
                case('f'): {
                    double dVal= va_arg( list, double );
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << fixed << dVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << fixed << dVal << " ";
                }
                break;
                case('e'): {
                    double dVal= va_arg( list, double );
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << scientific << dVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << scientific << dVal << " ";
                }
                break;
                case('s'): {
                    string sVal= va_arg( list, const char* );
                    if(LoggerDevice::Display == mvOutputDevice[mode])
                        cout << sVal << " ";
                    if(LoggerDevice::Disk == mvOutputDevice[mode])
                        ofs  << sVal << " ";
                }
                break;
                default:
                    break;
                }
            }
        };
        va_end( list );

        if(LoggerDevice::Display == mvOutputDevice[mode])
            cout << endl;
        if(LoggerDevice::Disk == mvOutputDevice[mode])
            ofs << endl;

    }//if ( nCurrentState ) end
}

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
