/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Logger.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
    mvOutputDevice[LoggerMode::Debug] = LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Error] = LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Warn]  = LoggerDevice::Display;
    mvOutputDevice[LoggerMode::Info]  = LoggerDevice::Display;
    mcFillStr = '_';
    msLoggerPref  = "Logger          ";
    mdDummyValue= 0.0;
    muDummyValue= 0;
    miDummyValue= 0;
}
CLogger::~CLogger()
{
}
void CLogger::InfoDisplay()
{
    cout << Info::Header()  << endl;
}
void CLogger::setFileBaseName(const string& filebase)
{
    msLogFileBaseName = filebase;
}
void CLogger::initializeLogFile(const uint& rank)
{
    myRank= rank;
    msLogFileName = msPathName + msLogFileBaseName + "_";
    msLogFileName += boost::lexical_cast<string>(myRank);
    msLogFileName += ".log";
    msOutputString = "Remove old LogFile";
    if(remove(msLogFileName.c_str()) == 0){
        cout << msLoggerPref << setfill(mcFillStr) << setw(60) << msOutputString << endl;
    }
    ofs.open(msLogFileName.c_str(), ios::out | ios::app);
}
void CLogger::finalizeLogFile()
{
    ofs.close();
}
void CLogger::setMode(const uint& mode)
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
            cout << msLoggerPref << setfill(mcFillStr) << setw(60) << "Logger Mode mismatch, at setMode()" << endl;
            break;
    }
}
void CLogger::setProperty(const uint& mode, const uint& outputDevice)
{
    if(outputDevice > LoggerDevice::Display){
        cout << msLoggerPref << setfill(mcFillStr) << setw(60) << "Logger Device mismatch, at setProperty()" << endl;
        return;
    }
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
            cout << msLoggerPref << setfill(mcFillStr) << setw(60) << "Logger Mode mismatch, at setProperty()" << endl;
            break;
    }    
}
bool CLogger::BooleanMode(const uint& mode)
{
    bool bCheck(false);
    if(mode >= LoggerMode::Invalid){
        cout << msLoggerPref << setfill(mcFillStr) << setw(60) << "Logger Mode mismatch, ";
        bCheck=false;
    }else{
        bCheck=true;
    }
    return bCheck;
}
string& CLogger::ModeToString(const uint& mode)
{
    switch(mode){
        case(LoggerMode::MWDebug):
            msModeString = "MW_Debug        ";
            break;
        case(LoggerMode::Debug):
            msModeString = "Debug           ";
            break;
        case(LoggerMode::Error):
            msModeString = "Error           ";
            break;
        case(LoggerMode::Warn):
            msModeString = "Warn            ";
            break;
        case(LoggerMode::Info):
            msModeString = "Info            ";
            break;
        default:
            msModeString = "LogMode Mismatch";
            break;
    }
    return msModeString;
}
void CLogger::Monitor(const uint& mode, const uint& id, const double& value, const string& message)
{
    if(!BooleanMode(mode)){
        cout <<  "at Monitor 1" << endl;
        return;
    }
    msOutputString = " ID=";
    ArgToString(id,value);
    msOutputString += " @" + message;
    if(mode >= mnCurrentState){
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
    }
}
void CLogger::Monitor(const uint& mode, const uint& id, const vdouble& vValue, const string& message)
{
    if(!BooleanMode(mode)){
        cout <<  "at Monitor 2" << endl;
        return;
    }
    msOutputString = " ID=";
    ArgToString(id,vValue);
    msOutputString += " @" + message;
    if(mode >= mnCurrentState){
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
    }
}
void CLogger::Monitor(const uint& mode, const uint& id, const vint& vValue, const string& message)
{
    if(!BooleanMode(mode)){
        cout <<  "at Monitor 3" << endl;
        return;
    }
    msOutputString = " ID=";
    ArgToString(id,vValue);
    msOutputString += " @" + message;
    if(mode >= mnCurrentState){
        uint i;
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
    }
}
string& CLogger::ArgToString(const uint& id, const double& value)
{
    msOutputString += boost::lexical_cast<string>(id);
    msOutputString += ",";
    msOutputString += boost::lexical_cast<string>(value) + ",";
    return msOutputString;
}
string& CLogger::ArgToString(const uint& id, const vdouble& vValue)
{
    uint i;
    msOutputString += boost::lexical_cast<string>(id);
    msOutputString += ",";
    for(i=0; i < vValue.size(); i++){
        msOutputString += boost::lexical_cast<string>(vValue[i]) + ",";
    };
    return msOutputString;
}
string& CLogger::ArgToString(const uint& id, const vint& vValue)
{
    uint i;
    msOutputString += boost::lexical_cast<string>(id);
    msOutputString += ",";
    for(i=0; i < vValue.size(); i++){
        msOutputString += boost::lexical_cast<string>(vValue[i]) + ",";
    };
    return msOutputString;
}
string& CLogger::InfoToString(const uint& num)
{
    msOutputString += boost::lexical_cast<string>(num);
    return msOutputString;
}
string& CLogger::InfoToString(const double& val)
{
    msOutputString += boost::lexical_cast<string>(val);
    return msOutputString;
}
void CLogger::Info(const uint& mode, const string& message)
{
    if(!BooleanMode(mode)){
        cout <<  "at Info()" << endl;
        return;
    }
    msOutputString = message;
    if(mode >= mnCurrentState){
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
    }
}
void CLogger::Info(const uint& mode, const string& message1, const string& message2)
{
    if(!BooleanMode(mode)){
        cout <<  "at Info()" << endl;
        return;
    }
    msOutputString = message1 + message2;
    if(mode >= mnCurrentState){
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
    }
}
void CLogger::Info(const uint& mode, const string& message, const uint& num)
{
    if(!BooleanMode(mode)){
        cout <<  "at Info()" << endl;
        return;
    }
    msOutputString = message + ", ";
    InfoToString(num);
    if(mode >= mnCurrentState){
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
    }
}
void CLogger::Info(const uint& mode, const string& message, const double& val)
{
    if(!BooleanMode(mode)){
        cout <<  "at Info()" << endl;
        return;
    }
    msOutputString = message + ", ";
    InfoToString(val);
    if(mode >= mnCurrentState){
        if(LoggerDevice::Display == mvOutputDevice[mode]){
            cout << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
        if(LoggerDevice::Disk == mvOutputDevice[mode]){
            ofs << ModeToString(mode) << setfill(mcFillStr) << setw(60) << msOutputString << endl;
        }
    }
}
double& CLogger::getDDummyValue()
{
    return mdDummyValue;
}
uint& CLogger::getUDummyValue()
{
    return muDummyValue;
}
int& CLogger::getIDummyValue()
{
    return miDummyValue;
}
