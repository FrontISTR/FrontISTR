//
//  GeneralBoundary.cpp
//
//
//
//                      2009.05.18
//                      2009.05.18
//                      k.Takeda
#include "GeneralBoundary.h"
#include "Logger.h"
using namespace pmw;

CGeneralBoundary::CGeneralBoundary()
{
    ;
}

CGeneralBoundary::~CGeneralBoundary()
{
//    //debug
//    cout << "~CGeneralBoundary" << endl;
}


// 境界値のセット for each val
//
void CGeneralBoundary::setValue(const double& val, const uint& i)
{
    if(i < mvValue.size()){
        mvValue[i]=val;
    }
}

void CGeneralBoundary::setValue(const vdouble& vVal)
{
    if(mvValue.size()==vVal.size()){
        mvValue = vVal;
    }else{
        Utility::CLogger *pLogger=Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"value.size() mismatch at Boundary::setValue");
    }   
}

// 境界値の提供 for each val
//
double& CGeneralBoundary::getValue(const uint& i)
{
    if(i < mvValue.size()){
        return mvValue[i];
    }else{
        Utility::CLogger *pLogger=Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"index mismatch at Boundary::getValue");
        
        return mvValue[0];
    }
}

