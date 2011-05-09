/* 
 * File:   DiscreteVolume.h
 * Author: ktakeda
 *
 * VolumeをTetraに分割するTree
 *
 *  スカラー三重積による体積計算用途
 *
 *  Prism-> 3個のTetra
 *  Hexa -> 2個のPrism -> 6個のTetra
 *
 * Created on 2010/06/22, 15:02
 */
#include "TypeDef.h"
#include "ElementType.h"

namespace pmw{
#ifndef _DISCRETEVOLUME_H
#define	_DISCRETEVOLUME_H
class CDiscreteVolume{
private:
    CDiscreteVolume();
public:
    static CDiscreteVolume* Instance(){
        static CDiscreteVolume moDisVol;
        return &moDisVol;
    }
    virtual ~CDiscreteVolume();

private:
    uint mPrismDiscre[3][4];// Prismの頂点番号で表した,3個のTetra
    uint mHexaDiscre[6][4]; // Hexaの頂点番号で表した,6個のTetra
public:
    uint* getPrismDiscrete(const uint& index){ return mPrismDiscre[index];}
    uint* getHexaDiscrete(const uint& index){ return mHexaDiscre[index];}
};
#endif	/* _DISCRETEVOLUME_H */
}

