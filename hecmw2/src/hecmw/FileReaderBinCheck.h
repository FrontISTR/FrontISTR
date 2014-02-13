/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBinCheck.h
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
#include "TypeDef.h"
#include "CommonFile.h"
#include "FileBlockName.h"
namespace FileIO
{
#ifndef FILEREADER_BIN_CHECK_H
#define	FILEREADER_BIN_CHECK_H
class CFileReaderBinCheck
{
private:
    CFileReaderBinCheck();
public:
    static CFileReaderBinCheck* Instance() {
        static CFileReaderBinCheck oReadHeader;
        return &oReadHeader;
    }
    virtual ~CFileReaderBinCheck();
private:
    bool mbLittleEndian;
    bool isLittleEndian_Sys();
    bool isBigEndian_Sys();
    static int32 mBit32;
    static int32 mBit64;
    bool mb32;
    bool mb64;
public:
    bool isLittleEndian_File() {
        return mbLittleEndian;
    }
    bool isBigEndian_File()   {
        return !mbLittleEndian;
    }
    bool isByteOrderSwap();
    void ByteOrderSwap(uiint& nVal);
    void ByteOrderSwap(iint& nVal);
    void ByteOrderSwap32(uint32& nVal);
    void ByteOrderSwap32(int32& nVal);
    void ByteOrderSwap64(uint64& nVal);
    void ByteOrderSwap64(int64& nVal);
    void ByteOrderSwap(double& dVal);
    bool Read_bin(ifstream& ifs);
    bool is32Bit() {
        return mb32;
    }
    bool is64Bit() {
        return mb64;
    }
};
#endif	/* FILEREADERHEADER_H */
}
