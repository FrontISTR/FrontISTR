//
//  Film.cpp
//
//
//
//                              2009.05.18
//                              2009.05.18
//                              k.Takeda

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
//    //debug
//    cout << "~CFilm" << endl;
}

// Accessor
//
void CFilm::setHeatTrans(const double& heat_trans)
{
    mHeatTrans = heat_trans;
}

double& CFilm::getHeatTrans()
{
    return mHeatTrans;
}

