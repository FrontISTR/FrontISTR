/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   SolidBase.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
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
