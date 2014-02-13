/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderCommNode_CM2.cpp
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
#include "FileReaderCommNode_CM2.h"
using namespace FileIO;
CFileReaderCommNodeCM2::CFileReaderCommNodeCM2()
{
    ;
}
CFileReaderCommNodeCM2::~CFileReaderCommNodeCM2()
{
    ;
}
string CFileReaderCommNodeCM2::Name()
{
    return  "FileReaderCommNodeCM2";
}


bool CFileReaderCommNodeCM2::Read(ifstream& ifs, string& sLine)
{
    uiint  mgLevel(0);
    uiint  comID, mesh_id, numOfCommNode;
    uiint  comm_node_id, node_id;
    vdouble vCoord;
    vCoord.resize(3);
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartCommNodeCM2()) ) {
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommNodeCM2", sLine);
        while(true) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommNodeCM2()) ) break;
            iss.clear();
            iss.str(sLine);
            iss >> comID >> mesh_id >> numOfCommNode;
            for(uiint i=0; i< numOfCommNode; i++) {
                sLine = getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> comm_node_id >> node_id >> vCoord[0] >> vCoord[1] >> vCoord[2];
                mpFactory->GeneCommNodeCM2(mgLevel, mesh_id, node_id, comID, comm_node_id, vCoord);
            };
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderCommNodeCM2::Read_bin(ifstream& ifs)
{
    return true;
}
