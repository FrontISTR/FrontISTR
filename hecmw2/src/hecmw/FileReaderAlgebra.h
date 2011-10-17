/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderAlgebra.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReader.h"
#include "FileReaderBinCheck.h" 
namespace FileIO{
#ifndef FILEREADERALGEBRA_H
#define	FILEREADERALGEBRA_H
class CFileReaderAlgebra:public CFileReader{
public:
    CFileReaderAlgebra();
    virtual ~CFileReaderAlgebra();
private:
    vuint mvAlgebraDOF;
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
    uiint getNumOfEquation(){ return mvAlgebraDOF.size();}
    uiint& getEquationDOF(const uiint& ieq){ return mvAlgebraDOF[ieq];}
};
#endif	/* FILEREADERALGEBRA_H */
}
