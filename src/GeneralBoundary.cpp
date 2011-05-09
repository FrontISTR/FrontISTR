/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   GeneralBoundary.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "GeneralBoundary.h"
#include "Logger.h"
using namespace pmw;
CGeneralBoundary::CGeneralBoundary()
{
    ;
}
CGeneralBoundary::~CGeneralBoundary()
{
}
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
