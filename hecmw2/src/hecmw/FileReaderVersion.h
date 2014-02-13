/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderVersion.h
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
#include "FileReader.h"

namespace FileIO
{
#ifndef VERSION_4685922d_8e72_42a5_a0b7_b5da00e61bd9
#define VERSION_4685922d_8e72_42a5_a0b7_b5da00e61bd9
class CFileReaderVersion:public CFileReader
{
public:
    CFileReaderVersion();
    virtual ~CFileReaderVersion();
private:
    string msVersion;
public:
    //現行バージョン:"04"
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
    virtual string Name();

    string& getVersion();
};
#endif //VERSION_4685922d_8e72_42a5_a0b7_b5da00e61bd9
}


