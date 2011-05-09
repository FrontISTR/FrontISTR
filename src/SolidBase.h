//
//
//
//				2008.11.27
//				2008.11.27
//				k.Takeda
#ifndef CSOLID_BASE_HH_50DA75CF_7692_46ce_9AEF_86A7C5877E70
#define CSOLID_BASE_HH_50DA75CF_7692_46ce_9AEF_86A7C5877E70

#include "Element.h"

namespace pmw{
class CSolidBase : public CElement
{
public:
    CSolidBase(void);
    virtual ~CSolidBase(void);

protected:
    static uint mnBaseType;

public:
    virtual const uint& getEntityType(){ return mnBaseType;}
};
}
#endif
