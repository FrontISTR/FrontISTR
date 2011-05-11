//
//
//              2008.12.1
//              2008.12.1
//              k.Takeda

#include "ShellBase.h"
#include "ElementType.h"
using namespace pmw;

uiint CShellBase::mnBaseType = BaseElementType::Shell;

//
//
CShellBase::CShellBase():CElement()
{
    ;
}

CShellBase::~CShellBase()
{
//    //debug
//    cout << "~CShellBase" << endl;
}
