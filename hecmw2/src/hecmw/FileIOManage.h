/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileIOConfig.h
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
#include "CommonFile.h"
#include "TypeDef.h"
#include "FileBlockName.h"
#include "Logger.h"

namespace FileIO
{
#ifndef Ee5a0f5d_4f8e_42ab_bea9_b9de5c26c31c
#define Ee5a0f5d_4f8e_42ab_bea9_b9de5c26c31c

class CFileIOManage
{
private:
    CFileIOManage();
public:
    virtual ~CFileIOManage();

    static CFileIOManage* Instance() {
        static CFileIOManage oManage;
        return &oManage;
    }

private:
    double mdVersion;
public:
    void setFileVersion(double version);
    double& getFileVersion();

    bool isAvailable_NumFrom();//境界条件数式が有効:書式Ver0.4以上
};

#endif
}

