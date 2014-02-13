/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/DiscreteVolume.h
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
#include "TypeDef.h"
#include "ElementType.h"
namespace pmw
{
#ifndef _DISCRETEVOLUME_H
#define	_DISCRETEVOLUME_H
class CDiscreteVolume
{
private:
    CDiscreteVolume();
public:
    static CDiscreteVolume* Instance() {
        static CDiscreteVolume moDisVol;
        return &moDisVol;
    }
    virtual ~CDiscreteVolume();
private:
    uiint mPrismDiscre[3][4];
    uiint mHexaDiscre[6][4];
public:
    uiint* getPrismDiscrete(const uiint& index) {
        return mPrismDiscre[index];
    }
    uiint* getHexaDiscrete(const uiint& index) {
        return mHexaDiscre[index];
    }
};
#endif	/* _DISCRETEVOLUME_H */
}
