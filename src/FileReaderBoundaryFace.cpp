//
//  FileReaderBoundaryFace.cpp
//
//
//
//                          2009.05.22
//                          2009.05.22
//                          k.Takeda
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

// method
// --
bool CFileReaderBoundaryFace::Read(ifstream& ifs, string& sLine)
{
    uint    nElemID, nFaceID, numOfBFace, nMeshID, maxID, minID, dof, bndType;
    uint    nMGLevel(0);
    vdouble vValue;
    string  s_bndType;

    if(TagCheck(sLine, FileBlockName::StartBoundaryFace()) ){

        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());

        // 境界Face数, MeshID, MaxID, MinID
        iss >> numOfBFace >> nMeshID >> maxID >> minID;
        //////////////// nMGLevel とりあえず"0"にしてある.
        mpFactory->reserveBoundaryFace(nMGLevel, nMeshID, numOfBFace);


        vValue.clear();
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryFace())) break;

            istringstream iss(sLine.c_str());
            // DOF, BoundaryType, ElementID,  FaceID
            iss >> dof >> s_bndType >> nElemID >> nFaceID;

            vValue.resize(dof);

            for(uint i=0; i< dof; i++) iss >> vValue[i];

            // ファイル中に記述されている境界タイプ文字列をunsigned int に変換
            if(s_bndType=="Pressure")       bndType = pmw::BoundaryTypeFace::Pressure;
            if(s_bndType=="TractionVector") bndType = pmw::BoundaryTypeFace::TractionVector;
            if(s_bndType=="Temp")         bndType = pmw::BoundaryTypeFace::Temp;
            if(s_bndType=="Thermal_Flux") bndType= pmw::BoundaryTypeFace::Thermal_Flux;

            // bndType の型判定によりBoundaryFaceにセットする値を判別
            mpFactory->GeneBoundaryFace(nMGLevel, nMeshID, nElemID, nFaceID, dof, bndType, vValue);

            mpLogger->Monitor(Utility::LoggerMode::MWDebug, nElemID, vValue, "FileReaderBoundaryFace");
        };
        return true;
    }else{
        return false;
    }
}


