/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Film.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Film.h"
#include "ElementType.h"
using namespace pmw;
uint CFilm::mnBaseType = BaseElementType::MaterialElement;
uint CFilm::mnElemType = ElementType::Polygon;
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
