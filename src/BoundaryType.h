/* 
 * File:   BoundaryType.h
 * Author: ktakeda
 *
 * Modify     2009/05/25
 * Created on 2009/05/14, 14:13
 */

#ifndef _BOUNDARYTYPE_H_f6a34053_61db_4c0c_8aeb_0a964eb6ee0f
#define	_BOUNDARYTYPE_H_f6a34053_61db_4c0c_8aeb_0a964eb6ee0f

namespace pmw{
union BoundaryTypeNode{
    enum {
        Load,
        Disp,
        Velo,
        Accel,
        Temp,
        Thermal_Flux
    };
};
union BoundaryTypeFace{
    enum {
        Pressure,
        TractionVector,
        Temp,
        Thermal_Flux
    };
};
union BoundaryTypeVolume{
    enum {
        Accel,
        Gravity,
        Centrifugal_Force,//遠心力
        Heat
    };
};

union BoundaryTypeFlow{
    enum {
        Velo_X,
        Velo_Y,
        Velo_Z,
        Pressure,
        Temperature,
        Heat,
        Thermal_Flux,
        OtherScalar
    };
};
}
#endif	/* _BOUNDARYTYPE_H */

