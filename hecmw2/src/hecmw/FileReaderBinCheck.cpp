//
//  FileReaderBinCheck.cpp
//
//
//              2011.03.18
//              k.Takeda
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


//bool CFileReaderBinCheck::Read(ifstream& ifs, string& sLine)
//{
//    Utility::CLogger *pLogger= Utility::CLogger::Instance();
//    pLogger->Info(Utility::LoggerMode::Warn, "invalid method, FileReaderBinCheck::Read(ifs, sLine) ");
//
//    return false;
//
//    istringstream iss;
//    string sToken;
//    short nVal;//2 byte
//    char *pCh; //1 byte
//
//    iss.clear();
//    iss.str(sLine);
//
//    while(iss){
//        iss >> sToken;
//
//        if(sToken=="@"){
//            cout << "Endian" << flush;//debug
//
//            iss >> nVal;//定数==1
//
//            pCh = (char*)&nVal;//先頭の1byte
//
//            if(*pCh){
//                mbLittleEndian = true;//Little endian
//                cout << " Little" << flush;//debug
//            }else{
//                mbLittleEndian = false;//Big endian
//                cout << " Big" << flush;//debug
//            }
//
//            uint nBit;
//            iss >> nBit;//ビット数
//
//            cout << " Bit" << nBit << endl;//debug
//
//            if( isByteOrderSwap() ) ByteOrderSwap(nBit);
//
//            return true;
//        }else{
//            return false;
//        }
//    };
//}

bool CFileReaderBinCheck::Read_bin(ifstream& ifs)
{
    uiint nLength= FileBlockName::BinCheck_Len();
    char cTag[nLength+1];
    
    ifs.read( cTag, nLength );

    cTag[nLength]='\0';
    string sTag=cTag;

    if( sTag == FileBlockName::StartBinCheck() ){
        
        int32 nNum;
        ifs.read( (char*)&nNum, sizeof(int32) );//定数=1 の読み込み

        char *pCh;
        pCh = (char*)&nNum;

        if(*pCh){
            mbLittleEndian = true; //リトル・エンディアン
        }else{
            mbLittleEndian = false;//ビッグ・エンディアン
        }
        
        int32 nBitNum;
        ifs.read((char*)&nBitNum, sizeof(int32));
        
        if(isByteOrderSwap()) ByteOrderSwap32(nBitNum);//バイト・オーダーが違えば反転.
        
        if(nBitNum==mBit32){ mb32=true; mb64=false;}
        if(nBitNum==mBit64){ mb64=true; mb32=false;}
        
        ifs.read(cTag, FileBlockName::End_Len());

        return true;
    }
    return false;
}

//--
// エンディアン判定  システム
//--
bool CFileReaderBinCheck::isLittleEndian_Sys() //コンピューターは リトル・エンディアン か
{
    short nVal=1;
    char *pCh;

    pCh = (char*)&nVal;

    if(*pCh){
        return true;//リトル・エンディアン
    }else{
        return false;//ビッグ・エンディアン
    }
}

bool CFileReaderBinCheck::isBigEndian_Sys() //コンピューターは ビッグ・エンディアン か
{
    short nVal=1;
    char *pCh;

    pCh = (char*)&nVal;

    if(*pCh){
        return false;//リトル・エンディアン
    }else{
        return true;//ビッグ・エンディアン
    }
}

// ----
// バイト・オーダーの入れ替え判定 : "入力ファイル"と"システム"のエンディアンの不一致なら"True"
// ----
bool CFileReaderBinCheck::isByteOrderSwap()
{
    if(isBigEndian_File() == isBigEndian_Sys()){
        return false;//バイト・オーダー そのまま
    }else{
        return  true;//バイト・オーダー 入れ替え
    }
}
// ----
// バイト・オーダーの入れ替え { バイト単位で、逆順に入れ替え }
// ----
void CFileReaderBinCheck::ByteOrderSwap(uiint& nVal)
{
    uiint dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );

    uiint nNumOfByte= sizeof(uiint);
    for(uiint iB=0; iB < nNumOfByte; iB++){
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap(iint& nVal)
{
    iint dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );

    uiint nNumOfByte= sizeof(iint);
    for(uiint iB=0; iB < nNumOfByte; iB++){
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap32(uint32& nVal)
{
    uint32 dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );

    uiint nNumOfByte= sizeof(uint32);
    for(uiint iB=0; iB < nNumOfByte; iB++){
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap32(int32& nVal)
{
    int32 dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );

    uiint nNumOfByte= sizeof(int32);
    for(uiint iB=0; iB < nNumOfByte; iB++){
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap64(uint64& nVal)
{
    uint64 dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );

    uiint nNumOfByte= sizeof(uint64);
    for(uiint iB=0; iB < nNumOfByte; iB++){
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap64(int64& nVal)
{
    int64 dTmp= nVal;
    char *pdTmp( (char*)&dTmp );
    char *pnVal( (char*)&nVal );

    uiint nNumOfByte= sizeof(int64);
    for(uiint iB=0; iB < nNumOfByte; iB++){
        pnVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}
void CFileReaderBinCheck::ByteOrderSwap(double& dVal)
{
    double dTmp= dVal;
    char *pdTmp( (char*)&dTmp );
    char *pdVal( (char*)&dVal );

    uiint nNumOfByte= sizeof(double);
    for(uiint iB=0; iB < nNumOfByte; iB++){
        pdVal[iB] = pdTmp[ nNumOfByte-1 -iB ];
    };
}




