/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBinCheck.cpp
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
#include "FileReaderBinCheck.h"
using namespace FileIO;
int32 CFileReaderBinCheck::mBit32 = 32;
int32 CFileReaderBinCheck::mBit64 = 64;
CFileReaderBinCheck::CFileReaderBinCheck()
{
    mbLittleEndian = false;
    mb32 = false;
    mb64 = false;
}
CFileReaderBinCheck::~CFileReaderBinCheck()
{
    ;
}
bool CFileReaderBinCheck::Read_bin(ifstream& ifs)
{
    uiint nLength= FileBlockName::BinCheck_Len();
    //char cTag[nLength+1];
    char* cTag;
    cTag = (char*)malloc(sizeof(char) * (nLength + 1));
    ifs.read( cTag, nLength );
    cTag[nLength]='\0';
    string sTag=cTag;
    if( sTag == FileBlockName::StartBinCheck() ) {
        int32 nNum;
        ifs.read( (char*)&nNum, sizeof(int32) );
        char *pCh;
        pCh = (char*)&nNum;
        if(*pCh) {
            mbLittleEndian = true;
        } else {
            mbLittleEndian = false;
        }
        int32 nBitNum;
        ifs.read((char*)&nBitNum, sizeof(int32));
        if(isByteOrderSwap()) ByteOrderSwap32(nBitNum);
        if(nBitNum==mBit32) {
            mb32=true;
            mb64=false;
        }
        if(nBitNum==mBit64) {
            mb64=true;
            mb32=false;
        }
        ifs.read(cTag, FileBlockName::End_Len());
        free(cTag);
        return true;
    }
    free(cTag);
    return false;
}
bool CFileReaderBinCheck::isLittleEndian_Sys()
{
    short nVal=1;
    char *pCh;
    pCh = (char*)&nVal;
    if(*pCh) {
        return true;
    } else {
        return false;
    }
}
bool CFileReaderBinCheck::isBigEndian_Sys()
{
    short nVal=1;
    char *pCh;
    pCh = (char*)&nVal;
    if(*pCh) {
        return false;
    } else {
        return true;
    }
}
bool CFileReaderBinCheck::isByteOrderSwap()
{
    if(isBigEndian_File() == isBigEndian_Sys()) {
        return false;
    } else {
        return  true;
    }
}
void CFileReaderBinCheck::ByteOrderSwap(uiint& nVal)
{
    uiint dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );
    uiint nNumOfByte= sizeof(uiint);
    for(uiint iB=0; iB < nNumOfByte; iB++) {
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap(iint& nVal)
{
    iint dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );
    uiint nNumOfByte= sizeof(iint);
    for(uiint iB=0; iB < nNumOfByte; iB++) {
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap32(uint32& nVal)
{
    uint32 dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );
    uiint nNumOfByte= sizeof(uint32);
    for(uiint iB=0; iB < nNumOfByte; iB++) {
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap32(int32& nVal)
{
    int32 dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );
    uiint nNumOfByte= sizeof(int32);
    for(uiint iB=0; iB < nNumOfByte; iB++) {
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap64(uint64& nVal)
{
    uint64 dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );
    uiint nNumOfByte= sizeof(uint64);
    for(uiint iB=0; iB < nNumOfByte; iB++) {
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap64(int64& nVal)
{
    int64 dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );
    uiint nNumOfByte= sizeof(int64);
    for(uiint iB=0; iB < nNumOfByte; iB++) {
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap(double& dVal)
{
    double dTmp= dVal;
    char *pdTmp( (char*)&dTmp );
    char *pdVal( (char*)&dVal );
    uiint nNumOfByte= sizeof(double);
    for(uiint iB=0; iB < nNumOfByte; iB++) {
        pdVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
