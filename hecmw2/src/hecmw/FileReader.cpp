/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReader.cpp
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
#include "FileReader.h"
using namespace FileIO;
CFileReader::CFileReader()
{
    mpLogger = Utility::CLogger::Instance();
    mpFileManage = CFileIOManage::Instance();
}
CFileReader::~CFileReader()
{
}
void CFileReader::setFactory(pmw::CMeshFactory *pFactory)
{
    mpFactory = pFactory;
}
bool CFileReader::TagCheck(string& s_line, const char* ctag)
{
    bool bcheck(false);
    uiint i;
    for(i=0; i< s_line.length(); i++) {
        if(s_line[i]=='\r') s_line[i]=' ';
        if(s_line[i]==',')  s_line[i]=' ';
    };
    istringstream iss(s_line.c_str());
    string stoken;
    while(iss >> stoken) {
        if(stoken==ctag) {
            bcheck=true;
            return bcheck;
        }
    };
    return bcheck;
}
string& CFileReader::getLineSt(ifstream& ifs)
{
    char c_Line[BUFFERLENGTH];
    ifs.getline(c_Line, sizeof(c_Line), '\n');
    msLine.clear();
    msLine = c_Line;
    uiint i;
    for(i=0; i< msLine.length(); i++) {
        if(msLine[i]=='\r') msLine[i]=' ';
        if(msLine[i]==',')  msLine[i]=' ';
        if(msLine[i]=='\t') msLine[i]=' ';
    };
    return msLine;
}
string& CFileReader::getLine(ifstream& ifs)
{
    msLine.clear();
    getline(ifs, msLine);
    uiint i;
    for(i=0; i< msLine.length(); i++) {
        if(msLine[i]=='\r') msLine[i]=' ';
        if(msLine[i]==',')  msLine[i]=' ';
        if(msLine[i]=='\t') msLine[i]=' ';
    };
    return msLine;
}
uiint CFileReader::IntElemType(string& sElemType)
{
    if(sElemType=="Hexa") {
        return pmw::ElementType::Hexa;
    } else if(sElemType=="Hexa2") {
        return pmw::ElementType::Hexa2;
    } else if(sElemType=="Tetra") {
        return pmw::ElementType::Tetra;
    } else if(sElemType=="Tetra2") {
        return pmw::ElementType::Tetra2;
    } else if(sElemType=="Prism") {
        return pmw::ElementType::Prism;
    } else if(sElemType=="Prism2") {
        return pmw::ElementType::Prism2;
    } else if(sElemType=="Quad") {
        return pmw::ElementType::Quad;
    } else if(sElemType=="Quad2") {
        return pmw::ElementType::Quad2;
    } else if(sElemType=="Triangle") {
        return pmw::ElementType::Triangle;
    } else if(sElemType=="Triangle2") {
        return pmw::ElementType::Triangle2;
    } else if(sElemType=="Beam") {
        return pmw::ElementType::Beam;
    } else if(sElemType=="Beam2") {
        return pmw::ElementType::Beam2;
    } else if(sElemType=="Point") {
        return pmw::ElementType::Point;
    } else {
        mpLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CFileReader::IntElemType");
        return mpLogger->getUDummyValue();
    }
}
uiint CFileReader::IntBndType(string& sBndType)
{
    if(sBndType=="Dirichlet") {
        return pmw::BoundaryType::Dirichlet;
    } else if(sBndType=="Neumann") {
        return pmw::BoundaryType::Neumann;
    } else {
        mpLogger->Info(Utility::LoggerMode::Error, sBndType, " invalid BoundaryType Name, CFileReader::IntBndType");
        return mpLogger->getUDummyValue();
    }
}
void CFileReader::Split(const string& s, char c, vstring& v)
{
    string::size_type i = 0;
    string::size_type j = s.find(c);
    while(j != string::npos) {
        v.push_back(s.substr(i, j-i));
        i = ++j;
        j = s.find(c, j);
        if(j == string::npos)
            v.push_back(s.substr(i, s.length()));
    };
}
string CFileReader::Cleaning(string& str)// 特殊記号の排除 => [0-9] | '.' | [A-Z] だけにする.
{
    string str2;
    for(string::size_type i=0; i < str.length(); i++) {
        if(isalnum(str[i]) || str[i]=='.') str2 += str[i];
    };

    return str2;
}
uiint CFileReader::getFileSize(ifstream& ifs)
{
    ifs.seekg(0, ios_base::end);
    uiint nSize = ifs.tellg();
    ifs.seekg(0, ios_base::beg);

    return nSize;
}
bool CFileReader::Check_End(ifstream& ifs)
{
    char cTok[4];
    ifs.read(cTok, 3);
    cTok[3]='\0';
    string sTok=cTok;
    if(sTok==FileBlockName::End()) {
        return true;
    } else {
        ifs.seekg(-3, ios_base::cur);
        return false;
    }
}
bool CFileReader::Check_IntSize(bool& b32, bool& bCheck, string& sClassName)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    b32= pBinCheck->is32Bit();
    bCheck=false;
    if(b32) {
        if(sizeof(uiint)==sizeof(uint32)) bCheck=true;
    } else {
        if(sizeof(uiint)==sizeof(uint64)) bCheck=true;
    }
    if(!bCheck) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        mpLogger->Info(Utility::LoggerMode::Error, sClassName+"::Read_bin, mismatch intger size");
        return false;
    }
    return true;
}
bool CFileReader::TagCheck_Bin(ifstream& ifs, bool bCheck, char cHead, const char *NameTag, const uiint& nLength)
{
    ifs.seekg(0, ios_base::beg);
    bCheck=false;
    string sTag;
    while(!ifs.eof()) {
        char cTok;
        ifs.read(&cTok, 1);
        if(cTok==cHead) {
            ifs.seekg(-1, ios_base::cur);
            //char cTag[nLength+1];
            char* cTag = new char[nLength+1];

            ifs.read(cTag, nLength);
            cTag[nLength]='\0';
            sTag=cTag;
            if( sTag==NameTag ) {
                bCheck=true;
                break;
            } else {
                ifs.seekg(-(nLength-1), ios_base::cur);
            }
            delete [] cTag;
        }
    };
    return bCheck;
}
short CFileReader::Read_ElementType(ifstream& ifs, char* cEName, short nLength, const char* cEType, string& sType)
{
    char cTail;
    ifs.read(cEName, nLength-1);
    cEName[nLength-1]='\0';
    sType=cEName;
    if(sType==cEType) {
        ifs.read(&cTail, 1);
        if(cTail=='2') {
            sType += "2";
            return 2;
        } else {
            ifs.seekg(-1, ios_base::cur);
            return 1;
        }
    }
    ifs.seekg(-nLength, ios_base::cur);
    sType="";
    return 0;
}
void CFileReader::Read_BndType(ifstream& ifs, string& sBndType)
{
    char cH;
    ifs.read(&cH, 1);
    if(cH=='D') {
        ifs.seekg(8, ios_base::cur);
        sBndType="Dirichlet";
    }
    if(cH=='N') {
        ifs.seekg(6, ios_base::cur);
        sBndType="Neumann";
    }
}
void CFileReader::Read_AnyName(ifstream& ifs, string& sName)
{
    sName.clear();
    int nOK;
    char cH;
    while(!ifs.eof()) {
        ifs.read(&cH, 1);
        nOK = isalnum(cH);
        if(nOK==0) {
            ifs.seekg(-1, ios_base::cur);
            break;
        } else {
            sName += cH;
        }
    };
}
