/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommFace.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderCommFace.h"
using namespace FileIO;
CFileReaderCommFace::CFileReaderCommFace()
{
    ;
}
CFileReaderCommFace::~CFileReaderCommFace()
{
    ;
}
bool CFileReaderCommFace::Read(ifstream& ifs, string& sLine)
{
    uint numOfCommFace;
    uint comID, meshID;
    uint nType, comm_face_id, elem_id, elem_ent_num;
    string sType;
    vuint vCommNodeID;
    uint mgLevel(0);
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartCommFace()) ){
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommFace", sLine);
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommFace()) ) break;
            iss.clear();
            iss.str(sLine);
            iss >> comID >> meshID >> numOfCommFace;
            uint iface;
            for(iface=0; iface< numOfCommFace; iface++){
                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> sType >> comm_face_id >> elem_id >> elem_ent_num;
                nType= IntElemType(sType);
                vCommNodeID.clear();
                switch(nType){
                    case(pmw::ElementType::Quad):
                        vCommNodeID.resize(4);
                        iss >> vCommNodeID[0] >> vCommNodeID[1] >> vCommNodeID[2] >> vCommNodeID[3];
                        break;
                    case(pmw::ElementType::Triangle):
                        vCommNodeID.resize(3);
                        iss >> vCommNodeID[0] >> vCommNodeID[1] >> vCommNodeID[2];
                        break;
                    case(pmw::ElementType::Beam):
                        vCommNodeID.resize(2);
                        iss >> vCommNodeID[0] >> vCommNodeID[1];
                        break;
                }
                mpFactory->GeneCommFace(mgLevel, comID, comm_face_id, meshID, elem_id, elem_ent_num, vCommNodeID);
            };
        };
        return true;
    }else{
        return false;
    }
}
