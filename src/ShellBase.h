/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ShellBase.h
|
|                     Written by T.Takeda,    2010/06/01
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
