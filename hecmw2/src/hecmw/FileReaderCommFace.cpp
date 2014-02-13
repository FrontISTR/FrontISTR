/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderCommFace.cpp
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
string CFileReaderCommFace::Name()
{
    return  "FileReaderCommFace";
}

bool CFileReaderCommFace::Read(ifstream& ifs, string& sLine)
{
    uiint numOfCommFace;
    uiint comID, meshID;
    uiint elem_type, comm_face_id, elem_id, elem_ent_num;
    string sType;
    vuint vCommNodeID;
    uiint mgLevel(0);
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartCommFace()) ) {
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommFace", sLine);
        while(true) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommFace()) ) break;
            iss.clear();
            iss.str(sLine);
            iss >> comID >> meshID >> numOfCommFace;
            uiint iface;
            for(iface=0; iface< numOfCommFace; iface++) {
                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                //形状, 通信面ID, 要素ID, 要素面番号(点の場合は局所節点番号)
                iss >> sType >> comm_face_id >> elem_id >> elem_ent_num;
                elem_type = IntElemType(sType);
                vCommNodeID.clear();
                switch(elem_type) {
                case(pmw::ElementType::Quad):
                    vCommNodeID.resize(4);
                    iss >> vCommNodeID[0] >> vCommNodeID[1] >> vCommNodeID[2] >> vCommNodeID[3];
                    break;
                case(pmw::ElementType::Quad2):
                    vCommNodeID.resize(8);
                    iss >> vCommNodeID[0] >> vCommNodeID[1] >> vCommNodeID[2] >> vCommNodeID[3]
                        >> vCommNodeID[4] >> vCommNodeID[5] >> vCommNodeID[6] >> vCommNodeID[7];
                    break;
                case(pmw::ElementType::Triangle):
                    vCommNodeID.resize(3);
                    iss >> vCommNodeID[0] >> vCommNodeID[1] >> vCommNodeID[2];
                    break;
                case(pmw::ElementType::Triangle2):
                    vCommNodeID.resize(6);
                    iss >> vCommNodeID[0] >> vCommNodeID[1] >> vCommNodeID[2]
                        >> vCommNodeID[3] >> vCommNodeID[4] >> vCommNodeID[5];
                    break;
                case(pmw::ElementType::Beam):
                    vCommNodeID.resize(2);
                    iss >> vCommNodeID[0] >> vCommNodeID[1];
                    break;
                case(pmw::ElementType::Beam2):
                    vCommNodeID.resize(3);
                    iss >> vCommNodeID[0] >> vCommNodeID[1] >> vCommNodeID[2];
                    break;
                case(pmw::ElementType::Point):
                    vCommNodeID.resize(1);
                    iss >> vCommNodeID[0];
                    break;
                }
                mpFactory->GeneCommFace(mgLevel, comID, comm_face_id, meshID, elem_id, elem_ent_num, elem_type, vCommNodeID);
            };
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderCommFace::Read_bin(ifstream& ifs)
{
    return true;
}
