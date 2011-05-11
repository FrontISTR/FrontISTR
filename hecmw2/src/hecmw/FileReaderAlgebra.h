/* 
 * File:   FileReaderAlgebra.h
 * Author: ktakeda
 *
 * Created on 2011/03/09, 14:34
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef FILEREADERALGEBRA_H
#define	FILEREADERALGEBRA_H
class CFileReaderAlgebra:public CFileReader{
public:
    CFileReaderAlgebra();
    virtual ~CFileReaderAlgebra();

private:
    vuint mvAlgebraDOF;//各 線形方程式のDOF

public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);

    uiint getNumOfEquation(){ return mvAlgebraDOF.size();}//線形方程式の個数
    uiint& getEquationDOF(const uiint& ieq){ return mvAlgebraDOF[ieq];}//各方程式のDOF
};
#endif	/* FILEREADERALGEBRA_H */
}


