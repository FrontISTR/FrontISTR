/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderAlgebra.h
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
#ifndef FILEREADERALGEBRA_H
#define	FILEREADERALGEBRA_H
class CFileReaderAlgebra:public CFileReader
{
public:
    CFileReaderAlgebra();
    virtual ~CFileReaderAlgebra();
private:
    uiint  mnNumOfLevel;
    vvuint mvAlgebraDOF;
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);

    virtual string Name();

    uiint getNumOfLevel() {
        return mnNumOfLevel;
    }
    uiint getNumOfEquation() {
        return mvAlgebraDOF.size();
    }
    uiint getNumOfParts() {
        return mvAlgebraDOF[0].size();
    }
    uiint& getEquationDOF(const uiint& ieq, const uiint& ipart) {
        return mvAlgebraDOF[ieq][ipart];
    }
};
#endif	/* FILEREADERALGEBRA_H */
}
