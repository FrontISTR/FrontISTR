/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderContactMesh.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderContactMesh.h"
using namespace FileIO;
CFileReaderContactMesh::CFileReaderContactMesh()
{
    ;
}
CFileReaderContactMesh::~CFileReaderContactMesh()
{
    ;
}
bool CFileReaderContactMesh::Read(ifstream& ifs, string& sLine)
{
    uint numOfContactMesh, maxID, minID;
    uint my_rank, trans_rank;
    uint contactID, meshID, nodeID, rank, elemID, elemFaceID, shapeType;
    bool bmesh(false);
    uint maslave;
    uint mgLevel(0);
    uint conNodeID, numOfConNode, contactFaceID, numOfContactFace;
    vuint vConNodeID;
    vuint quadConNodeID; quadConNodeID.resize(4);
    vuint triConNodeID;  triConNodeID.resize(3);
    vuint beamConNodeID; beamConNodeID.resize(2);
    vuint pointConNodeID;pointConNodeID.resize(1);
    vdouble vCoord; vCoord.resize(3);
    string shapeStr;
    vstring strData; strData.resize(3);
    string sParamType;
    uint numOfDisp,numOfScalar;
    if(TagCheck(sLine, FileBlockName::StartContactMesh()) ){
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndContactMesh()) ) break; 
            istringstream iss(sLine.c_str());
            iss >> numOfContactMesh >> maxID >> minID;
            uint icont;
            for(icont=0; icont< numOfContactMesh; icont++){
                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> contactID >> my_rank >> trans_rank;
                mpFactory->GeneContactMesh(contactID, my_rank, trans_rank);
                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> numOfConNode >> maxID >> minID;
                uint icnode;
                for(icnode=0; icnode< numOfConNode; icnode++){
                    sLine= getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);
                    iss >> conNodeID >> vCoord[0] >> vCoord[1] >> vCoord[2] >> strData[0] >> strData[1] >> rank >> maslave;
                    iss >> sParamType >> numOfDisp >> numOfScalar;
                    bmesh=false;
                    if(strData[0]!="-" && strData[1]!="-"){
                        bmesh= true;
                        meshID = boost::lexical_cast<unsigned int>(strData[0]);
                        nodeID = boost::lexical_cast<unsigned int>(strData[1]);
                    }
                    mpFactory->GeneContactNode(mgLevel, contactID, conNodeID, vCoord, 
                                                sParamType, numOfDisp, numOfScalar,
                                                bmesh, meshID, nodeID, rank, maslave);
                };
                for(maslave=0; maslave< 2; maslave++){
                    sLine= getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);
                    iss >> numOfContactFace >> maxID >> minID;
                    uint iface;
                    for(iface=0; iface< numOfContactFace; iface++){
                        sLine= getLineSt(ifs);
                        iss.clear();
                        iss.str(sLine);
                        iss >> contactFaceID >> strData[0]  >> strData[1] >> strData[2] >> shapeStr;
                        bmesh=false;
                        if(strData[0]!="-" && strData[1]!="-" && strData[2]!="-"){
                            bmesh= true;
                            meshID =  boost::lexical_cast<unsigned int>(strData[0]);
                            elemID =  boost::lexical_cast<unsigned int>(strData[1]);
                            elemFaceID = boost::lexical_cast<unsigned int>(strData[2]);
                        }
                        vConNodeID.clear();
                        if(shapeStr=="Quad"){
                            shapeType= pmw::ElementType::Quad;
                            iss >> quadConNodeID[0] >> quadConNodeID[1] >> quadConNodeID[2] >> quadConNodeID[3];
                            vConNodeID= quadConNodeID;
                        }
                        if(shapeStr=="Triangle"){
                            shapeType= pmw::ElementType::Triangle;
                            iss >> triConNodeID[0] >> triConNodeID[1] >> triConNodeID[2];
                            vConNodeID= triConNodeID;
                        }
                        if(shapeStr=="Beam"){
                            shapeType= pmw::ElementType::Beam;
                            iss >> beamConNodeID[0] >> beamConNodeID[1];
                            vConNodeID= beamConNodeID;
                        }
                        if(shapeStr=="Point"){
                            shapeType= pmw::ElementType::Point;
                            iss >> pointConNodeID[0];
                            vConNodeID= pointConNodeID;
                        }
                        uint face_rank;
                        iss >> face_rank;
                        switch(maslave){
                            case(0):
                                mpFactory->GeneMasterFace(contactID, shapeType, contactFaceID, bmesh, 
                                        meshID,elemID,elemFaceID,vConNodeID, face_rank);
                                break;
                            case(1):
                                mpFactory->GeneSlaveFace(contactID,shapeType,contactFaceID, bmesh, 
                                        meshID,elemID,elemFaceID, vConNodeID, face_rank);
                                break;
                        }
                    };
                };
            };
        };
        return true;
    }else{
        return false;
    }
}
