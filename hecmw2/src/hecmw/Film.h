/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Film.h
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
#include "CommonStd.h"
#include "TypeDef.h"
////#include "Element.h"

namespace pmw
{
#ifndef _FILM_H_ad201417_2f36_4da7_8847_40ffded43b0c
#define	_FILM_H_ad201417_2f36_4da7_8847_40ffded43b0c
class CFilm
{
public:
    CFilm();
    virtual ~CFilm();

private:
    vdouble mvTransCoeff;//---- 熱伝達係数 配列：線形方程式番号に対応

public:
    // 線形方程式別の伝達率
    void  initTransCoeff(const uiint& nNumOfEquation);
    uiint getNumOfEquation() {
        return mvTransCoeff.size();
    }

    void    setTransCoeff(const double& transCoeff, const uiint& ieq) {
        mvTransCoeff[ieq]= transCoeff;
    }
    double& getTransCoeff(const uiint& ieq) {
        return mvTransCoeff[ieq];
    }
};
#endif	/* _FILM_H */
}

