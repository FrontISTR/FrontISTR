//
//  FileReaderCommMesh.cpp
//
//
//
//                  2009.09.17
//                  2009.09.17
//                  k.Takeda
#include "FileReaderCommMesh.h"
using namespace FileIO;

// construct && destruct
//
CFileReaderCommMesh::CFileReaderCommMesh()
{
    ;
}

CFileReaderCommMesh::~CFileReaderCommMesh()
{
    ;
}

// CommMesh(通信領域の読み込み)
//
bool CFileReaderCommMesh::Read(ifstream& ifs, string& sLine)
{
    uint mgLevel(0);// mgLevel=0 ::ファイル入力時のマルチグリッド・レベルは恒等的に==0
    
    uint  numOfMesh, maxMeshID, minMeshID;
    uint  nMeshID, numOfCommMesh;
    uint  nCommMeshID, myRank,nTransmitRank;

    // Logger表示用
    string   sLinePara;
    vstring svLinePara; svLinePara.resize(3);
    string   white(" ");

    // MeshごとのCommMeshよみこみ
    // --
    if(TagCheck(sLine, FileBlockName::StartCommMesh()) ){
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommMesh", sLine);

        // CommMesh数,MaxID,MinID
        //
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfMesh >> maxMeshID >> minMeshID;

        // CommMeshオブジェクトの配列確保
        // --
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommMesh()) ) break;

            istringstream iss(sLine.c_str());

            // MeshID, CommMesh数 (MeshごとのCommMesh数)
            // --
            iss >> nMeshID >> numOfCommMesh;

            // CommMeshの配列確保
            // --
            mpFactory->reserveCommMesh(mgLevel, nMeshID, numOfCommMesh);//ファイル読み込みなので,mgLevel==0

            // CommMesh数ぶんのCommMeshID, myRank, TrasmiRank
            // --
            for(uint i=0; i< numOfCommMesh; i++){
                sLine = getLineSt(ifs);
                istringstream issin(sLine.c_str());

                issin >> nCommMeshID >> myRank >> nTransmitRank;
                
                // Loggerによる表示のための処理
                svLinePara[0]= boost::lexical_cast<string>(nCommMeshID);
                svLinePara[1]= boost::lexical_cast<string>(myRank);
                svLinePara[2]= boost::lexical_cast<string>(nTransmitRank);
                sLinePara = svLinePara[0] + white + svLinePara[1] + white + svLinePara[2];
                // --
                mpLogger->Info(Utility::LoggerMode::MWDebug, "CommMesh, myRank, TransmitRank => ", sLinePara);

                // 通信領域番号とランク(自分のRank,送信先Rank)
                mpFactory->GeneCommMesh(mgLevel, nMeshID, nCommMeshID, myRank, nTransmitRank);
            };
        };
        return true;
    }else{
        return false;
    }
}



