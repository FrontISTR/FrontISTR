/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderAlgebra.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
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
    bool b32, bCheck;
    string sClassName("FileReaderAlgebra");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;
    uiint nLength= FileBlockName::Algebra_Len();
    //char cTag[nLength+1];
	char* cTag = (char*)malloc(sizeof(char) * (nLength + 1));
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
        ifs.read(cTag, FileBlockName::End_Len());
		free(cTag);
        return true;
    }
	free(cTag);
    return false;
}
