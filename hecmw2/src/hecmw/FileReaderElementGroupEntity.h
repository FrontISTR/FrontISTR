/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderElementGroupEntity.h
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
#include "FileReaderBinCheck.h"
namespace FileIO
{
#ifndef FILEREADER_ELEMENTGROUP_ENTITY_H
#define	FILEREADER_ELEMENTGROUP_ENTITY_H
class CFileReaderElementGroupEntity:public CFileReader
{
public:
    CFileReaderElementGroupEntity();
    virtual ~CFileReaderElementGroupEntity();
public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);

    virtual string Name();
};
#endif	/* FILEREADERELEMENTGROUPENTITY_H */
}
