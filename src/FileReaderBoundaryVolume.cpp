/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderBoundaryVolume.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderBoundaryVolume.h"
using namespace FileIO;
CFileReaderBoundaryVolume::CFileReaderBoundaryVolume()
{
    ;
}
CFileReaderBoundaryVolume::~CFileReaderBoundaryVolume()
{
    ;
}
bool CFileReaderBoundaryVolume::Read(ifstream& ifs, string& sLine)
{
    uint    nElemID, numOfBVol, nMeshID, maxID, minID, dof, bndType;
    uint    nMGLevel(0);
    vdouble vValue;
    string  s_bndType;
    if(TagCheck(sLine, FileBlockName::StartBoundaryVolume()) ){
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfBVol >> nMeshID >> maxID >> minID;
        mpFactory->reserveBoundaryVolume(nMGLevel, nMeshID, numOfBVol);
        vValue.clear();
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryVolume())) break;
            istringstream iss(sLine.c_str());
            iss >> dof >> s_bndType >> nElemID;
            vValue.resize(dof);
            for(uint i=0; i< dof; i++) iss >> vValue[i];
            if(s_bndType=="Accel")            bndType = pmw::BoundaryTypeVolume::Accel;
            if(s_bndType=="Gravity")          bndType = pmw::BoundaryTypeVolume::Gravity;
            if(s_bndType=="Centrifugal_Force") bndType = pmw::BoundaryTypeVolume::Centrifugal_Force;
            if(s_bndType=="Heat")              bndType= pmw::BoundaryTypeVolume::Heat;
            mpFactory->GeneBoundaryVolume(nMGLevel, nMeshID, nElemID, dof, bndType, vValue);
            mpLogger->Monitor(Utility::LoggerMode::MWDebug, nElemID, vValue, "FileReaderBoundaryNode");
        };
        return true;
    }else{
        return false;
    }
}
