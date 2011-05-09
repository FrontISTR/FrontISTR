/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BeamBase.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef BEAM_BASE_HH_79610CD9_A0C2_4446_82F0_BFC537CF5373
#define BEAM_BASE_HH_79610CD9_A0C2_4446_82F0_BFC537CF5373
#include "CommonStd.h"
#include "Element.h"
namespace pmw{
class CBeamBase:public CElement
{
public:
    CBeamBase(void);
    virtual ~CBeamBase(void);
protected:
    static uint mnBaseType;
public:
    virtual const uint& getEntityType(){ return mnBaseType;}
};
}
#endif
