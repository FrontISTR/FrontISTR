//
//  FileReaderCommNode.cpp
//
//
//
//                  2009.09.18
//                  2009.09.18
//                  k.Takeda
#include "FileReaderCommNode.h"
using namespace FileIO;


// construct & destruct
//
CFileReaderCommNode::CFileReaderCommNode()
{
    ;
}
CFileReaderCommNode::~CFileReaderCommNode()
{
    ;
}

// method
//
bool CFileReaderCommNode::Read(ifstream& ifs, string& sLine)
{
    uint numOfNode, nCommMeshID, maxID, minID;
    uint nCommNodeID, nMeshID, nNodeID, nRank;

    uint  mgLevel(0); // MultiGrid Level==0 ::ファイル入力時は 0

    // CommMesh内のmvNode => MeshのNodeを取得してセット
    if(TagCheck(sLine, FileBlockName::StartCommNode()) ){
        //mpLogger->Info(Utility::LoggerMode::MWDebug, sLine);

        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());

        // CommMesh内ノード数, MeshID, CommMeshID, MaxID, MinID
        iss >> numOfNode >> nMeshID >> nCommMeshID >> maxID >> minID;

        // mvNodeの配列数確保(CommNodeというものは,存在しない.)
        // --
        mpFactory->reserveCommNode(mgLevel,nMeshID, nCommMeshID, numOfNode);

        // CommNode 読み込み => Mesh::Nodeからオブジェクト取得
        //   => CommMesh::mvNodeにセットされる.
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommNode()) ) break;

            istringstream iss(sLine.c_str());
            // CommNodeID, NodeID, rank
            iss >> nCommNodeID >> nNodeID >> nRank;

            // MeshからNodeを取得(CommNodeというものは,存在しない.),rankによりSend,Recv判定
            // --
            mpFactory->GeneCommNode(mgLevel, nCommNodeID, nMeshID, nCommMeshID, nNodeID, nRank);
        };
        return true;
    }else{
        return false;
    }
}


















