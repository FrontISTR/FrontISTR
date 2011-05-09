/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommNode_CM2.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
bool CFileReaderCommNodeCM2::Read(ifstream& ifs, string& sLine)
{
    uint  mgLevel(0);
    uint  comID, mesh_id, numOfCommNode;
    uint  comm_node_id, node_id;
    vdouble vCoord;  vCoord.resize(3);
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartCommNodeCM2()) ){
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommNodeCM2", sLine);
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommNodeCM2()) ) break;
            iss.clear();
            iss.str(sLine);
            iss >> comID >> mesh_id >> numOfCommNode;
            for(uint i=0; i< numOfCommNode; i++){
                sLine = getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> comm_node_id >> node_id >> vCoord[0] >> vCoord[1] >> vCoord[2];
                mpFactory->GeneCommNodeCM2(mgLevel, mesh_id, node_id, comID, comm_node_id, vCoord);
            };
        };
        return true;
    }else{
        return false;
    }
}
