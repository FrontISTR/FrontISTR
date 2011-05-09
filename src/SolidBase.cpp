//
//
//                          2008.11.27
//                          k.Takeda
#include "SolidBase.h"
#include "ElementType.h"
using namespace pmw;

uint CSolidBase::mnBaseType = BaseElementType::Solid;

CSolidBase::CSolidBase(void):CElement()
{
}

CSolidBase::~CSolidBase(void)
{
//    //debug
//    cout << "~CSolidBase" << endl;
}
