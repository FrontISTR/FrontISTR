/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReader.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
    uint i;
    for(i=0; i< s_line.length(); i++){
        if(s_line[i]=='\r') s_line[i]=' ';
        if(s_line[i]==',')  s_line[i]=' ';
    };
    istringstream iss(s_line.c_str());
    string stoken;
    while(iss >> stoken){
        if(stoken==ctag){
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
    uint i;
    for(i=0; i< msLine.length(); i++){
        if(msLine[i]=='\r') msLine[i]=' ';
        if(msLine[i]==',')  msLine[i]=' ';
    };
    return msLine;
}
uint CFileReader::IntElemType(string& sElemType)
{
    if(sElemType=="Hexa") return pmw::ElementType::Hexa;
    if(sElemType=="Tetra")return pmw::ElementType::Tetra;
    if(sElemType=="Prism")return pmw::ElementType::Prism;
    if(sElemType=="Pyramid")return pmw::ElementType::Pyramid;
    if(sElemType=="Quad") return pmw::ElementType::Quad;
    if(sElemType=="Triangle")return pmw::ElementType::Triangle;
    if(sElemType=="Beam") return pmw::ElementType::Beam;
}
