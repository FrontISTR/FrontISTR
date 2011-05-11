//
// FileReaderAlgebra.cpp
//
//              2011.03.09
//              k.Takeda
#include "FileReaderAlgebra.h"
using namespace FileIO;


CFileReaderAlgebra::CFileReaderAlgebra()
{
    ;
}
CFileReaderAlgebra::~CFileReaderAlgebra()
{
    ;
}

bool CFileReaderAlgebra::Read(ifstream& ifs, string& sLine)
{
    istringstream iss;

    // Algebra ブロック (リスタート・データ)
    if(TagCheck(sLine, FileBlockName::StartAlgebra()) ){
        sLine = getLineSt(ifs);

        iss.clear();
        iss.str(sLine);

        uiint nNumOfEquation;
        iss >> nNumOfEquation;

        
        mvAlgebraDOF.resize(nNumOfEquation);

        for(uiint ieq=0; ieq < nNumOfEquation; ieq++){
            sLine = getLineSt(ifs);
            iss.clear();
            iss.str(sLine);

            uiint nEqID;
            iss >> nEqID >> mvAlgebraDOF[ieq];
        };

        return true;
    }else{
        return false;
    }
}

bool CFileReaderAlgebra::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    
    //BinCheckのサイズ指定との整合性
    bool b32, bCheck;
    string sClassName("FileReaderAlgebra");

    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;

    
    uiint nLength= FileBlockName::Algebra_Len();
    char cTag[nLength+1];

    ifs.read(cTag, nLength);
    cTag[nLength]='\0';

    string sTag=cTag;

    if( sTag==FileBlockName::StartAlgebra() ){
        
        uiint nNumOfEquation;
        ifs.read((char*)&nNumOfEquation, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfEquation);

        mvAlgebraDOF.resize(nNumOfEquation);

        for(uiint ieq=0; ieq < nNumOfEquation; ieq++){
            uiint nEqID;
            ifs.read((char*)&nEqID, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nEqID);

            uiint nDOF;
            ifs.read((char*)&nDOF, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nDOF);
            mvAlgebraDOF[ieq]=nDOF;
        };
        ifs.read(cTag, FileBlockName::End_Len());//End

        return true;
    }
    return false;
}








