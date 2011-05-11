/* 
 * File:   FileWriterRes.h
 * Author: ktakeda
 *
 * Created on 2011/03/04, 16:58
 */
#include "FileWriter.h"

namespace pmw{
class CAssyModel;
class CMesh;
class CNode;
}

namespace FileIO{
#ifndef FILEWRITERRES_H
#define	FILEWRITERRES_H
class CFileWriterRes:public CFileWriter{
public:
    CFileWriterRes();
    virtual ~CFileWriterRes();

public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);//Non use.

    void WriteRes(ofstream& ofs);    //Resブロック出力(SolutionVectorの値)
    void WriteAlgebra(ofstream& ofs);//Algebraブロック出力(各線形方程式のDOF)

    void WriteRes_bin(ofstream& ofs);
    void WriteAlgebra_bin(ofstream& ofs);
};
#endif	/* FILEWRITERRES_H */
}

