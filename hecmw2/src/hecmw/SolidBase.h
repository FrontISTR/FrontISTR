/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/SolidBase.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
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
namespace pmw
{
class CSolidBase : public CElement
{
public:
    CSolidBase(void);
    virtual ~CSolidBase(void);
protected:
    static uiint mnBaseType;
public:
    virtual const uiint& getEntityType() {
        return mnBaseType;
    }
};
}
#endif
