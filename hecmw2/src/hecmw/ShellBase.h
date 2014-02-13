/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShellBase.h
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
#ifndef SHELL_BASE_HH_F8AEA60_C7E1_4b86_B9D0_23BB7B56185E
#define SHELL_BASE_HH_F8AEA60_C7E1_4b86_B9D0_23BB7B56185E
#include "CommonStd.h"
#include "Element.h"
namespace pmw
{
class CShellBase:public CElement
{
public:
    CShellBase(void);
    virtual ~CShellBase(void);
protected:
    static uiint mnBaseType;
public:
    virtual const uiint& getEntityType() {
        return mnBaseType;
    }
};
}
#endif
