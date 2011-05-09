/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   GeneralBoundary.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _GENERALBOUNDARY_H_efcafcfa_20ee_47d5_a9ef_cef945a71900
#define	_GENERALBOUNDARY_H_efcafcfa_20ee_47d5_a9ef_cef945a71900
#include "CommonStd.h"
#include "TypeDef.h"
#include "BoundaryType.h"
namespace pmw{
class CGeneralBoundary{
public:
    CGeneralBoundary();
    virtual ~CGeneralBoundary();
protected:
    uint mnType;
    vdouble mvValue;
public:
    virtual void initialize(const uint& bndType, const uint& dof)=0;
    const uint& getType(){ return mnType;}
    void setValue(const double& val, const uint& i);
    void setValue(const vdouble& vVal);
    double& getValue(const uint& i);
    vdouble& getValue(){ return mvValue;}
    uint numOfValue(){ return mvValue.size();}
};
}
#endif	/* _GENERALBOUNDARY_H */
