/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterVTK.h
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
#include "AssyVector.h"
namespace FileIO
{
#ifndef FILEWRITERVTK_H
#define	FILEWRITERVTK_H
class CFileWriterVTK:public CFileWriter
{
private:
    CFileWriterVTK();
public:
    static CFileWriterVTK* Instance() {
        static CFileWriterVTK moWriterVTK;
        return &moWriterVTK;
    }
    virtual ~CFileWriterVTK();
private:
    vector<vstring> mvLabel;
    vector<map<string, string> >  mmUnit;
    vector<map<string, uiint> >   mmDOF;
    vector<map<string, vdouble> > mmVecVal;
    vuint mvNumOfParameter;

    void WriteMesh(ofstream& ofs, const uiint& nNumParam, const uiint& iMesh, const uiint& mgLevel);
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
    void recVTK_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recVTK_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    void WriteFEM_Rank0(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel, string& basename, const uiint& nNumOfProcs);
    void WriteFEM(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel);
};
#endif	/* FILEWRITERVTK_H */
}
