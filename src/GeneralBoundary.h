/* 
 * File:   GeneralBoundary.h
 * Author: ktakeda
 *
 * Created on 2009/05/18, 18:15
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
    // Boundary Type 別の初期化
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

