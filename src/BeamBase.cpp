/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BeamBase.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "BeamBase.h"
#include "ElementType.h"
using namespace pmw;
uint CBeamBase::mnBaseType = BaseElementType::Beam;
CBeamBase::CBeamBase():CElement()
{
}
CBeamBase::~CBeamBase()
{
}
