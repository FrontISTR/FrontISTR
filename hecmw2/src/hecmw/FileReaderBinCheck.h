/* 
 * File:   FileReaderBinCheck.h
 * Author: ktakeda
 *
 * Created on 2011/03/18, 2:42
 */
#include "TypeDef.h"
#include "CommonFile.h"
#include "FileBlockName.h"


namespace FileIO{
#ifndef FILEREADER_BIN_CHECK_H
#define	FILEREADER_BIN_CHECK_H
class CFileReaderBinCheck{
private:
    CFileReaderBinCheck();
public:
    static CFileReaderBinCheck* Instance(){
        static CFileReaderBinCheck oReadHeader;
        return &oReadHeader;
    }
    virtual ~CFileReaderBinCheck();
    
private:
    bool mbLittleEndian;

    bool isLittleEndian_Sys();//コンピューターは リトル・エンディアン か
    bool isBigEndian_Sys();   //コンピューターは ビッグ・エンディアン か

    static int32 mBit32;
    static int32 mBit64;

    bool mb32;
    bool mb64;

public:
    bool isLittleEndian_File(){ return mbLittleEndian;} //入力ファイルは リトル・エンディアン か
    bool isBigEndian_File()   { return !mbLittleEndian;}//入力ファイルは ビッグ・エンディアン か

    bool isByteOrderSwap();//バイト・オーダーの入れ替え判定:"入力ファイル"と"システム"のエンディアンの不一致なら"True"

    // バイト・オーダーの入れ替え
    //
    void ByteOrderSwap(uiint& nVal);
    void ByteOrderSwap(iint& nVal);
    void ByteOrderSwap32(uint32& nVal);
    void ByteOrderSwap32(int32& nVal);
    void ByteOrderSwap64(uint64& nVal);
    void ByteOrderSwap64(int64& nVal);
    void ByteOrderSwap(double& dVal);

    // バイナリ読み込み
    bool Read_bin(ifstream& ifs);

    // 整数( iint,uiint )のビット数
    //
    bool is32Bit(){ return mb32;}
    bool is64Bit(){ return mb64;}
};
#endif	/* FILEREADERHEADER_H */
}


