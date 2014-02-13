/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileIOManage.cpp
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

#include "FileIOManage.h"
using namespace FileIO;

CFileIOManage::CFileIOManage()
{
    mdVersion = 0.0;
}
CFileIOManage::~CFileIOManage()
{
}
void CFileIOManage::setFileVersion(double version)
{
    mdVersion= version;
}
double& CFileIOManage::getFileVersion()
{
    return mdVersion;
}
//--
// 境界条件数式が有効:書式Ver0.4以上
//--
bool CFileIOManage::isAvailable_NumFrom()
{
    if(mdVersion >= 0.4) {
        return true;
    } else {
        return false;
    }
}

