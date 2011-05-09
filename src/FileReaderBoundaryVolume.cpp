//
//  FileReaderBoundaryVolume.cpp
//
//
//
//                          2009.05.22
//                          2009.05.22
//                          k.Takeda
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


// method
// --
bool CFileReaderBoundaryVolume::Read(ifstream& ifs, string& sLine)
{
    uint    nElemID, numOfBVol, nMeshID, maxID, minID, dof, bndType;
    uint    nMGLevel(0);
    vdouble vValue;
    string  s_bndType;

    //
    //
    if(TagCheck(sLine, FileBlockName::StartBoundaryVolume()) ){

        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());

        // 境界Volume数, MeshID, MaxID, MinID
        iss >> numOfBVol >> nMeshID >> maxID >> minID;
        //////////////// nMGLevel とりあえず"0"にしてある.
        mpFactory->reserveBoundaryVolume(nMGLevel, nMeshID, numOfBVol);


        vValue.clear();
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryVolume())) break;

            istringstream iss(sLine.c_str());
            // DOF, BoundaryType, ElementID
            iss >> dof >> s_bndType >> nElemID;

            vValue.resize(dof);

            for(uint i=0; i< dof; i++) iss >> vValue[i];

            // ファイル中に記述されている境界タイプ文字列をunsigned int に変換
            if(s_bndType=="Accel")            bndType = pmw::BoundaryTypeVolume::Accel;
            if(s_bndType=="Gravity")          bndType = pmw::BoundaryTypeVolume::Gravity;
            if(s_bndType=="Centrifugal_Force") bndType = pmw::BoundaryTypeVolume::Centrifugal_Force;//遠心力
            if(s_bndType=="Heat")              bndType= pmw::BoundaryTypeVolume::Heat;

            // bndType の型判定によりBoundaryVolumeにセットする値を判別
            mpFactory->GeneBoundaryVolume(nMGLevel, nMeshID, nElemID, dof, bndType, vValue);

            mpLogger->Monitor(Utility::LoggerMode::MWDebug, nElemID, vValue, "FileReaderBoundaryNode");
        };
        return true;
    }else{
        return false;
    }
}


