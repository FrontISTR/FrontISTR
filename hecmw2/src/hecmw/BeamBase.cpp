/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BeamBase.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include "BeamBase.h"
#include "ElementType.h"
using namespace pmw;
uiint CBeamBase::mnBaseType = BaseElementType::Beam;
CBeamBase::CBeamBase():CElement()
{
}
CBeamBase::~CBeamBase()
{
}
