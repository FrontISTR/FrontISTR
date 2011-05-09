//
//  CBeamBase CElement
//
//						2008.12.1
//						2008.12.1
//						k.Takeda
#ifndef BEAM_BASE_HH_79610CD9_A0C2_4446_82F0_BFC537CF5373
#define BEAM_BASE_HH_79610CD9_A0C2_4446_82F0_BFC537CF5373


#include "CommonStd.h"
#include "Element.h"

namespace pmw{
class CBeamBase:public CElement
{
public:
    CBeamBase(void);
    virtual ~CBeamBase(void);

protected:
    static uint mnBaseType;

public:
    virtual const uint& getEntityType(){ return mnBaseType;}
};
}
#endif
