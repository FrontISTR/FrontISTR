/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Film.cpp
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
#include "HEC_MPI.h"
#include "Film.h"
#include "ElementType.h"
using namespace pmw;
uiint CFilm::mnBaseType = BaseElementType::MaterialElement;
uiint CFilm::mnElemType = ElementType::Polygon;
CFilm::CFilm()
{
    ;
}
CFilm::~CFilm()
{
}
void CFilm::setHeatTrans(const double& heat_trans)
{
    mHeatTrans = heat_trans;
}
double& CFilm::getHeatTrans()
{
    return mHeatTrans;
}
