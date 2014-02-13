/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderAlgebra.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderAlgebra.h"
using namespace FileIO;
CFileReaderAlgebra::CFileReaderAlgebra()
{
    mnNumOfLevel= 0;
}
CFileReaderAlgebra::~CFileReaderAlgebra()
{
    ;
}
string CFileReaderAlgebra::Name()
{
    return "FileReaderAlgebra";
}
bool CFileReaderAlgebra::Read(ifstream& ifs, string& sLine)
{
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartAlgebra()) ) {
        sLine = getLineSt(ifs);
        iss.clear();
        iss.str(sLine);

        // Level数、線形方程式数、Mesh数
        uiint nNumOfEquation, nNumOfMesh;
        iss >> mnNumOfLevel >> nNumOfEquation >> nNumOfMesh;

        ////cout << "Level数:" << mnNumOfLevel << " 方程式数:" << nNumOfEquation << " Mesh数:" << nNumOfMesh << endl;//debug

        mvAlgebraDOF.resize(nNumOfEquation);
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++) {
            mvAlgebraDOF[ieq].resize(nNumOfMesh);
        };

        // 方程式番号、Mesh番号、DOF
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++) {
            sLine = getLineSt(ifs);
            iss.clear();
            iss.str(sLine);
            uiint nEq;
            iss >> nEq;
            if(nEq != ieq) mpLogger->Info(Utility::LoggerMode::Error, "FileReaderAlgebra::Read, nEq != ieq");

            for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
                uiint nMesh, nDOF;
                iss >> nMesh >> nDOF;

                ////cout << " Mesh:" << nMesh << " DOF:" << nDOF << endl;//debug

                if(nMesh != imesh) mpLogger->Info(Utility::LoggerMode::Error, "FileReaderAlgebra::Read, nPart != ipart");

                mvAlgebraDOF[ieq][imesh]= nDOF;
            };
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderAlgebra::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();

    bool bOrder= pBinCheck->isByteOrderSwap();
    bool b32, bCheck;
    string sClassName("FileReaderAlgebra");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;

    uiint nLength= FileBlockName::Algebra_Len();

    //char cTag[nLength+1];
    char* cTag = (char*)malloc(sizeof(char) * (nLength + 1));
    ifs.read(cTag, nLength);

    cTag[nLength]='\0';
    string sTag=cTag;

    if( sTag==FileBlockName::StartAlgebra() ) {
        uiint nNumOfEquation, nNumOfMesh;

        // Level数、線形方程式数、Mesh数
        ifs.read((char*)&mnNumOfLevel, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(mnNumOfLevel);
        ifs.read((char*)&nNumOfEquation, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(nNumOfEquation);
        ifs.read((char*)&nNumOfMesh, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(nNumOfMesh);

        ////cout << "Binary: Level数:" << mnNumOfLevel << " 方程式数:" << nNumOfEquation << " Mesh数:" << nNumOfMesh << endl;//debug

        mvAlgebraDOF.resize(nNumOfEquation);
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++) {
            mvAlgebraDOF[ieq].resize(nNumOfMesh);
        };

        // 方程式番号、Mesh番号、DOF
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++) {
            uiint nEq;
            ifs.read((char*)&nEq, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(nEq);

            for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
                uiint nMesh, nDOF;
                ifs.read((char*)&nMesh, sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(nMesh);
                ifs.read((char*)&nDOF, sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(nDOF);

                ////cout << "Binary: Mesh:" << nMesh << " DOF:" << nDOF << endl;//debug

                mvAlgebraDOF[ieq][imesh]=nDOF;
            };
        };
        ifs.read(cTag, FileBlockName::End_Len());
        free(cTag);
        return true;
    }
    free(cTag);
    return false;
}
