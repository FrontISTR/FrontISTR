/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderContactMesh.cpp
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
#include "FileReaderContactMesh.h"
using namespace FileIO;
using namespace boost;
CFileReaderContactMesh::CFileReaderContactMesh()
{
    ;
}
CFileReaderContactMesh::~CFileReaderContactMesh()
{
    ;
}
string CFileReaderContactMesh::Name()
{
    return  "FileReaderContactMesh";
}

bool CFileReaderContactMesh::Read(ifstream& ifs, string& sLine)
{
    uiint nNumOfContactMesh, maxID, minID;
    uiint my_rank;//------------- myRank
    uiint nProp(0);//------------ 接合面属性

    uiint contactID, meshID, nodeID, rank, elemID, elemFaceID, shapeType;

    bool  bmesh(false);
    uiint maslave;
    uiint mgLevel(0);

    uiint conNodeID, numOfConNode, contactFaceID, numOfContactFace;
    vuint vConNodeID;
    vuint quadConNodeID;
    quadConNodeID.resize(4);
    vuint quad2ConNodeID;
    quad2ConNodeID.resize(8);
    vuint triConNodeID;
    triConNodeID.resize(3);
    vuint tri2ConNodeID;
    tri2ConNodeID.resize(6);
    vuint beamConNodeID;
    beamConNodeID.resize(2);
    vuint beam2ConNodeID;
    beam2ConNodeID.resize(3);
    vuint pointConNodeID;
    pointConNodeID.resize(1);

    vdouble vCoord;
    vCoord.resize(3);
    string shapeStr;
    vstring strData;
    strData.resize(3);
    string sParamType;
    uiint numOfDisp,numOfScalar;

    if(TagCheck(sLine, FileBlockName::StartContactMesh()) ) {
        while(true) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndContactMesh()) ) break;
            istringstream iss(sLine.c_str());
            iss >> nNumOfContactMesh >> maxID >> minID;
            uiint icont;
            for(icont=0; icont< nNumOfContactMesh; icont++) {
                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);

                //--
                // 接合メッシュID,  myRank, 属性番号
                //--
                char_separator<char> sep(" \t\n");
                tokenizer< char_separator<char> > tokens(sLine, sep);
                uiint nCount(0);
                typedef tokenizer< char_separator<char> >::iterator Iter;
                for(Iter it=tokens.begin(); it != tokens.end(); ++it) {
                    string str = *it;
                    if(nCount==0) {
                        contactID  = atoi(str.c_str());   //
                    }
                    if(nCount==1) {
                        my_rank    = atoi(str.c_str());   //
                    }
                    if(nCount==2) {
                        nProp      = atoi(str.c_str());   //
                    }

                    nCount++;
                };

                mpFactory->GeneContactMesh(contactID, my_rank, nProp);

                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> numOfConNode >> maxID >> minID;
                uiint icnode;
                //--
                // ContactNode
                //--
                for(icnode=0; icnode< numOfConNode; icnode++) {
                    sLine= getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);
                    iss >> conNodeID >> vCoord[0] >> vCoord[1] >> vCoord[2] >> strData[0] >> strData[1] >> rank >> maslave;
                    iss >> sParamType >> numOfDisp >> numOfScalar;
                    bmesh=false;
                    if(strData[0]!="-" && strData[1]!="-") {
                        bmesh= true;
////////////////////    meshID = boost::lexical_cast<unsigned int>(strData[0]);
////////////////////    nodeID = boost::lexical_cast<unsigned int>(strData[1]);
                        istringstream sstoken;
                        sstoken.str(strData[0]);
                        sstoken >> meshID;
                        sstoken.clear();
                        sstoken.str(strData[1]);
                        sstoken >> nodeID;
                    }
                    mpFactory->GeneContactNode(mgLevel, contactID, conNodeID, vCoord,
                                               sParamType, numOfDisp, numOfScalar,
                                               bmesh, meshID, nodeID, rank, maslave);
                };
                //--
                // ContactFace
                //--
                for(maslave=0; maslave< 2; maslave++) {
                    sLine= getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);
                    iss >> numOfContactFace >> maxID >> minID;
                    uiint iface;
                    for(iface=0; iface< numOfContactFace; iface++) {
                        sLine= getLineSt(ifs);
                        iss.clear();
                        iss.str(sLine);
                        iss >> contactFaceID >> strData[0]  >> strData[1] >> strData[2] >> shapeStr;
                        bmesh=false;
                        if(strData[0]!="-" && strData[1]!="-" && strData[2]!="-") {
                            bmesh= true;
////////////////////////    meshID =  boost::lexical_cast<unsigned int>(strData[0]);
////////////////////////    elemID =  boost::lexical_cast<unsigned int>(strData[1]);
////////////////////////    elemFaceID = boost::lexical_cast<unsigned int>(strData[2]);
                            istringstream sstoken;
                            sstoken.str(strData[0]);
                            sstoken >> meshID;
                            sstoken.clear();
                            sstoken.str(strData[1]);
                            sstoken >> elemID;
                            sstoken.clear();
                            sstoken.str(strData[2]);
                            sstoken >> elemFaceID;
                        }
                        vConNodeID.clear();
                        shapeType = IntElemType(shapeStr);
                        switch(shapeType) {
                        case(pmw::ElementType::Quad):
                            iss >> quadConNodeID[0] >> quadConNodeID[1] >> quadConNodeID[2] >> quadConNodeID[3];
                            vConNodeID= quadConNodeID;
                            break;
                        case(pmw::ElementType::Quad2):
                            iss >> quad2ConNodeID[0] >> quad2ConNodeID[1] >> quad2ConNodeID[2] >> quad2ConNodeID[3]
                                >> quad2ConNodeID[4] >> quad2ConNodeID[5] >> quad2ConNodeID[6] >> quad2ConNodeID[7];
                            vConNodeID= quad2ConNodeID;
                            break;
                        case(pmw::ElementType::Triangle):
                            iss >> triConNodeID[0] >> triConNodeID[1] >> triConNodeID[2];
                            vConNodeID= triConNodeID;
                            break;
                        case(pmw::ElementType::Triangle2):
                            iss >> tri2ConNodeID[0] >> tri2ConNodeID[1] >> tri2ConNodeID[2]
                                >> tri2ConNodeID[3] >> tri2ConNodeID[4] >> tri2ConNodeID[5];
                            vConNodeID= tri2ConNodeID;
                            break;
                        case(pmw::ElementType::Beam):
                            iss >> beamConNodeID[0] >> beamConNodeID[1];
                            vConNodeID= beamConNodeID;
                            break;
                        case(pmw::ElementType::Beam2):
                            iss >> beam2ConNodeID[0] >> beam2ConNodeID[1] >> beam2ConNodeID[2];
                            vConNodeID= beam2ConNodeID;
                            break;
                        case(pmw::ElementType::Point):
                            iss >> pointConNodeID[0];
                            vConNodeID= pointConNodeID;
                            break;
                        default:
                            mpLogger->Info(Utility::LoggerMode::Error, "Not exist ShapeType, FileReaderContactMesh::Read");
                            break;
                        }
                        uiint face_rank;
                        iss >> face_rank;
                        switch(maslave) {
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
    } else {
        return false;
    }
}
bool CFileReaderContactMesh::Read_bin(ifstream& ifs)
{
    return true;
}
