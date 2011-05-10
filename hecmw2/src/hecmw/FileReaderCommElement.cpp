//
//  FileReaderCommElement.cpp
//
//
//
//                  2009.09.18
//                  2009.09.18
//                  k.Takeda
#include "FileReaderCommElement.h"
using namespace FileIO;

// construct & destruct
//
CFileReaderCommElement::CFileReaderCommElement()
{
    ;
}
CFileReaderCommElement::~CFileReaderCommElement()
{
    ;
}

// method
//
bool CFileReaderCommElement::Read(ifstream& ifs, string& sLine)
{
    uint mgLevel(0);

    uint  numOfCommElement, nMeshID, nCommMeshID, nMaxCommID, nMinCommID;
    uint  nElementID;
    vuint vCommNodeID;

    string sElemType;
    uint   nElemType;

    // CommElement
    if(TagCheck(sLine, FileBlockName::StartCommElement()) ){
        //mpLogger->Info(Utility::LoggerMode::MWDebug, sLine);

        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());

        // CommElement数, MeshID, CommMeshID, MaxID, MinID
        iss >> numOfCommElement >> nMeshID >> nCommMeshID >> nMaxCommID >> nMinCommID;
        
        
        // CommElement配列確保
        // --
        mpFactory->reserveCommElement(mgLevel, nMeshID, nCommMeshID, numOfCommElement);

        // CommElement 読み込み => Mesh::Elementを取得してCommElementにセット
        //
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommElement()) ) break;
            
            istringstream iss(sLine.c_str());
            // 書式：要素タイプ, ElementID, 各頂点のCommNodeID(下のswitch文で読み込み)
            iss >> sElemType >> nElementID;

            nElemType= IntElemType(sElemType);//文字列をunsigned int に変換

            //各頂点のCommNodeID(CommMeshのノード・インデックス)
            vCommNodeID.clear();
            uint ivert;
            switch(nElemType){
                case(pmw::ElementType::Hexa):
                    vCommNodeID.resize(8);
                    for(ivert=0; ivert< 8; ivert++) iss >> vCommNodeID[ivert];
                    break;
                case(pmw::ElementType::Tetra):
                    vCommNodeID.resize(4);
                    for(ivert=0; ivert< 4; ivert++) iss >> vCommNodeID[ivert];
                    break;
                case(pmw::ElementType::Prism):
                    vCommNodeID.resize(6);
                    for(ivert=0; ivert< 6; ivert++) iss >> vCommNodeID[ivert];
                    break;
//                case(pmw::ElementType::Pyramid):
//                    vCommNodeID.resize(5);
//                    for(ivert=0; ivert< 5; ivert++) iss >> vCommNodeID[ivert];
//                    break;
                case(pmw::ElementType::Quad):
                    vCommNodeID.resize(4);
                    for(ivert=0; ivert< 4; ivert++) iss >> vCommNodeID[ivert];
                    break;
                case(pmw::ElementType::Triangle):
                    vCommNodeID.resize(3);
                    for(ivert=0; ivert< 3; ivert++) iss >> vCommNodeID[ivert];
                    break;
                case(pmw::ElementType::Beam):
                    vCommNodeID.resize(2);
                    for(ivert=0; ivert< 2; ivert++) iss >> vCommNodeID[ivert];
                    break;
                default:
                    break;
            }
            
            // CommElementを生成,Elementをセット,要素の頂点のランクをセット(Factory内で,頂点のCommNodeIDからランクを取得している)
            // --
            mpFactory->GeneCommElement(mgLevel, nMeshID, nCommMeshID, nElemType, nElementID, vCommNodeID);
        };
        return true;
    }else{
        return false;
    }
}










