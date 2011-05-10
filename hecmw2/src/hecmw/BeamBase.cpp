//
//
//			2008.12.1
//			2008.12.1
//			k.Takeda

#include "BeamBase.h"
#include "ElementType.h"
using namespace pmw;

uint CBeamBase::mnBaseType = BaseElementType::Beam;
//
//
CBeamBase::CBeamBase():CElement()
{
}

CBeamBase::~CBeamBase()
{
    ////debug
    //cout << "~CBeamBase" << endl;
}

