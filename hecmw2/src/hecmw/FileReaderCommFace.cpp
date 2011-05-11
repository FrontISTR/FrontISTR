//
//  FileReaderCommFace.cpp
//
//
//
//
//              2010.03.12
//              k.Takeda
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
    uiint numOfCommFace;
    uiint comID, meshID;
    uiint elem_type, comm_face_id, elem_id, elem_ent_num;
    string sType;
    vuint vCommNodeID;
    
    uiint mgLevel(0);
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartCommFace()) ){
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommFace", sLine);
        
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommFace()) ) break;
            iss.clear();
            iss.str(sLine);

            iss >> comID >> meshID >> numOfCommFace;

            uiint iface;
            for(iface=0; iface< numOfCommFace; iface++){
                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                
                iss >> sType >> comm_face_id >> elem_id >> elem_ent_num;

                elem_type = IntElemType(sType);// 文字列”Quad,...”等をElementTypeに変換
                
                vCommNodeID.clear();
                switch(elem_type){
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
                }
                
                mpFactory->GeneCommFace(mgLevel, comID, comm_face_id, meshID, elem_id, elem_ent_num, elem_type, vCommNodeID);
            };
        };
        return true;
    }else{
        return false;
    }
}


bool CFileReaderCommFace::Read_bin(ifstream& ifs)
{
    return true;
}













