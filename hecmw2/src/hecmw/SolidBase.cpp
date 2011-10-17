/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/SolidBase.cpp
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
#include "SolidBase.h"
#include "ElementType.h"
using namespace pmw;
uiint CSolidBase::mnBaseType = BaseElementType::Solid;
CSolidBase::CSolidBase(void):CElement()
{
}
CSolidBase::~CSolidBase(void)
{
}
