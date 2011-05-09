/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BoundaryType.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _BOUNDARYTYPE_H_f6a34053_61db_4c0c_8aeb_0a964eb6ee0f
#define	_BOUNDARYTYPE_H_f6a34053_61db_4c0c_8aeb_0a964eb6ee0f
namespace pmw{
struct BoundaryTypeNode{
    enum {
        Load,
        Disp,
        Velo,
        Accel,
        Temp,
        Thermal_Flux
    };
};
struct BoundaryTypeFace{
    enum {
        Pressure,
        TractionVector,
        Temp,
        Thermal_Flux
    };
};
struct BoundaryTypeVolume{
    enum {
        Accel,
        Gravity,
        Centrifugal_Force,
        Heat
    };
};
struct BoundaryTypeFlow{
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
