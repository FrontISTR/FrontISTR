/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterRes.h
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
#include "FileWriter.h"
namespace pmw
{
class CAssyModel;
class CMesh;
class CNode;
}
namespace FileIO
{
#ifndef FILEWRITERRES_H
#define	FILEWRITERRES_H
class CFileWriterRes:public CFileWriter
{
public:
    CFileWriterRes();
    virtual ~CFileWriterRes();
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
    void WriteRes(ofstream& ofs);
    void WriteAlgebra(ofstream& ofs);
    void WriteRes_bin(ofstream& ofs);
    void WriteAlgebra_bin(ofstream& ofs);
};
#endif	/* FILEWRITERRES_H */
}
