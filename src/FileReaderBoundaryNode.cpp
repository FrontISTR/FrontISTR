//
//  FileReaderBoundaryNode.cpp
//
//
//
//                          2009.05.22
//                          2009.05.22
//                          k.Takeda
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

    // 
    //
    if(TagCheck(sLine, FileBlockName::StartBoundaryNode()) ){

        //mpLogger->Info(Utility::LoggerMode::MWDebug, "BoundaryNode", sLine);//debug

        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());

        // 境界Node数, MeshID, MaxID, MinID, DOF, BoundaryID
        iss >> numOfBNode >> nMeshID >> maxID >> minID;
        //////////////// nMGLevel とりあえず"0"にしてある.
        mpFactory->reserveBoundaryNode(nMGLevel, nMeshID, numOfBNode);


        vValue.clear();
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryNode())) break;

            istringstream iss(sLine.c_str());
            // DOF, BoundaryType, NodeID
            iss >> dof >> s_bndType >> nNodeID;

            vValue.resize(dof);

            for(uint i=0; i< dof; i++) iss >> vValue[i];

            // ファイル中に記述されている境界タイプ文字列をunsigned int に変換
            if(s_bndType=="Load") bndType = pmw::BoundaryTypeNode::Load;
            if(s_bndType=="Disp") bndType = pmw::BoundaryTypeNode::Disp;
            if(s_bndType=="Velo") bndType = pmw::BoundaryTypeNode::Velo;
            if(s_bndType=="Accel")        bndType= pmw::BoundaryTypeNode::Accel;
            if(s_bndType=="Temp")         bndType = pmw::BoundaryTypeNode::Temp;
            if(s_bndType=="Thermal_Flux") bndType= pmw::BoundaryTypeNode::Thermal_Flux;

            // bndType の型判定によりBoundaryNodeにセットする値を判別
            mpFactory->GeneBoundaryNode(nMGLevel, nMeshID, nNodeID, dof, bndType, vValue);

            mpLogger->Monitor(Utility::LoggerMode::MWDebug, nNodeID, vValue, "FileReaderBoundaryNode");
        };
        return true;
    }else{
        return false;
    }
}
















