//
//
//
//						2008.12.1
//						k.Takeda
#ifndef SHELL_BASE_HH_F8AEA60_C7E1_4b86_B9D0_23BB7B56185E
#define SHELL_BASE_HH_F8AEA60_C7E1_4b86_B9D0_23BB7B56185E

#include "CommonStd.h"
#include "Element.h"

namespace pmw{
class CShellBase:public CElement{
public:
    CShellBase(void);
    virtual ~CShellBase(void);

protected:
    static uint mnBaseType;

public:
    virtual const uint& getEntityType(){ return mnBaseType;}
};
}
#endif
