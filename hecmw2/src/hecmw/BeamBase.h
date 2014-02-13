/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BeamBase.h
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
#ifndef BEAM_BASE_HH_79610CD9_A0C2_4446_82F0_BFC537CF5373
#define BEAM_BASE_HH_79610CD9_A0C2_4446_82F0_BFC537CF5373
#include "CommonStd.h"
#include "Element.h"
namespace pmw
{
class CBeamBase:public CElement
{
public:
    CBeamBase(void);
    virtual ~CBeamBase(void);
protected:
    static uiint mnBaseType;
public:
    virtual const uiint& getEntityType() {
        return mnBaseType;
    }
};
}
#endif
