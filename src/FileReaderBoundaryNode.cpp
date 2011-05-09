/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderBoundaryNode.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderBoundaryNode.h"
using namespace FileIO;
CFileReaderBoundaryNode::CFileReaderBoundaryNode()
{
    ;
}
CFileReaderBoundaryNode::~CFileReaderBoundaryNode()
{
    ;
}
bool CFileReaderBoundaryNode::Read(ifstream& ifs, string& sLine)
{
    uint    nNodeID, numOfBNode, nMeshID, maxID, minID, dof, bndType;
    uint    nMGLevel(0);
    vdouble vValue;
    string  s_bndType;
    if(TagCheck(sLine, FileBlockName::StartBoundaryNode()) ){
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfBNode >> nMeshID >> maxID >> minID;
        mpFactory->reserveBoundaryNode(nMGLevel, nMeshID, numOfBNode);
        vValue.clear();
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryNode())) break;
            istringstream iss(sLine.c_str());
            iss >> dof >> s_bndType >> nNodeID;
            vValue.resize(dof);
            for(uint i=0; i< dof; i++) iss >> vValue[i];
            if(s_bndType=="Load") bndType = pmw::BoundaryTypeNode::Load;
            if(s_bndType=="Disp") bndType = pmw::BoundaryTypeNode::Disp;
            if(s_bndType=="Velo") bndType = pmw::BoundaryTypeNode::Velo;
            if(s_bndType=="Accel")        bndType= pmw::BoundaryTypeNode::Accel;
            if(s_bndType=="Temp")         bndType = pmw::BoundaryTypeNode::Temp;
            if(s_bndType=="Thermal_Flux") bndType= pmw::BoundaryTypeNode::Thermal_Flux;
            mpFactory->GeneBoundaryNode(nMGLevel, nMeshID, nNodeID, dof, bndType, vValue);
            mpLogger->Monitor(Utility::LoggerMode::MWDebug, nNodeID, vValue, "FileReaderBoundaryNode");
        };
        return true;
    }else{
        return false;
    }
}
