//
//  FileReaderCommMesh2.cpp
//
//
//
//
//          2010.03.12
//          k.Takeda
#include "FileReaderCommMesh2.h"
using namespace FileIO;

CFileReaderCommMesh2::CFileReaderCommMesh2()
{
    ;
}
CFileReaderCommMesh2::~CFileReaderCommMesh2()
{
    ;
}

bool CFileReaderCommMesh2::Read(ifstream& ifs, string& sLine)
{
    uint  mgLevel(0);// mgLevel=0 ::ファイル入力時のマルチグリッド・レベルは恒等的に==0
    uint  nMeshID, numOfCommMesh;
    uint  nCommMeshID, numOfFace, numOfCommNode, myRank, nTransmitRank;
    
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartCommMesh2()) ){
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommMesh2", sLine);
        
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommMesh2()) ) break;
            iss.clear();
            iss.str(sLine);

            iss >> nMeshID >> numOfCommMesh;// MeshIDごとの"CommMesh2"数 (Meshごとの"CommMesh2"数)

            //// CommMeshの配列確保
            //// --
            //mpFactory->reserveCommMesh2(mgLevel, nMeshID, numOfCommMesh);//ファイル読み込みなので,mgLevel==0

            // "CommMesh2"数ぶんのCommMeshID, myRank, TrasmiRank
            for(uint i=0; i< numOfCommMesh; i++){
                sLine = getLineSt(ifs);
                iss.clear();
                iss.str(sLine);

                iss >> nCommMeshID >> numOfFace >> numOfCommNode >> myRank >> nTransmitRank;
                
                // 通信領域番号とランク(自分のRank,送信先Rank)
                mpFactory->GeneCommMesh2(mgLevel, nMeshID, nCommMeshID,
                                            numOfFace, numOfCommNode, myRank, nTransmitRank);
            };
        };
        return true;
    }else{
        return false;
    }
}













