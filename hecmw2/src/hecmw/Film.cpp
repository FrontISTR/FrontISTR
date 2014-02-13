/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Film.cpp
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
#include "Film.h"

using namespace pmw;
////uiint CFilm::mnBaseType = BaseElementType::MaterialElement;
////uiint CFilm::mnElemType = ElementType::Polygon;
CFilm::CFilm()
{
    ;
}
CFilm::~CFilm()
{
}
//--
// initialize 熱伝達率
//--
void CFilm::initTransCoeff(const uiint& nNumOfEquation)
{
    mvTransCoeff.resize(nNumOfEquation);

    // デフォルト:伝達率= 1.0
    for(uiint ieq=0; ieq < nNumOfEquation; ieq++)
        mvTransCoeff[ieq]= 1.000;
}



