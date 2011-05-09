/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderBoundaryFace.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderBoundaryFace.h"
using namespace FileIO;
CFileReaderBoundaryFace::CFileReaderBoundaryFace()
{
    ;
}
CFileReaderBoundaryFace::~CFileReaderBoundaryFace()
{
    ;
}
bool CFileReaderBoundaryFace::Read(ifstream& ifs, string& sLine)
{
    uint    nElemID, nFaceID, numOfBFace, nMeshID, maxID, minID, dof, bndType;
    uint    nMGLevel(0);
    vdouble vValue;
    string  s_bndType;
    if(TagCheck(sLine, FileBlockName::StartBoundaryFace()) ){
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfBFace >> nMeshID >> maxID >> minID;
        mpFactory->reserveBoundaryFace(nMGLevel, nMeshID, numOfBFace);
        vValue.clear();
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryFace())) break;
            istringstream iss(sLine.c_str());
            iss >> dof >> s_bndType >> nElemID >> nFaceID;
            vValue.resize(dof);
            for(uint i=0; i< dof; i++) iss >> vValue[i];
            if(s_bndType=="Pressure")       bndType = pmw::BoundaryTypeFace::Pressure;
            if(s_bndType=="TractionVector") bndType = pmw::BoundaryTypeFace::TractionVector;
            if(s_bndType=="Temp")         bndType = pmw::BoundaryTypeFace::Temp;
            if(s_bndType=="Thermal_Flux") bndType= pmw::BoundaryTypeFace::Thermal_Flux;
            mpFactory->GeneBoundaryFace(nMGLevel, nMeshID, nElemID, nFaceID, dof, bndType, vValue);
            mpLogger->Monitor(Utility::LoggerMode::MWDebug, nElemID, vValue, "FileReaderBoundaryFace");
        };
        return true;
    }else{
        return false;
    }
}
