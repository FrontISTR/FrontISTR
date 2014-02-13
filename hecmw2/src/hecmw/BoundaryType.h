/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryType.h
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
namespace pmw
{
#ifndef _BOUNDARYTYPE_H_
#define	_BOUNDARYTYPE_H_
struct BoundaryType {
    enum {
        Dirichlet,
        Neumann,
        Limit
    };
};
struct BoundaryMeshType {
    enum {
        Node,
        Edge,
        Face,
        Volume,
        Limit
    };
};
struct BoundaryTypeNode {
    enum {
        Load,
        Disp,
        Velo,
        Accel,
        Temp,
        Thermal_Flux,
        Limit
    };
};
struct BoundaryTypeFace {
    enum {
        Pressure,
        TractionVector,
        Temp,
        Thermal_Flux,
        Limit
    };
};
struct BoundaryTypeVolume {
    enum {
        Accel,
        Gravity,
        Centrifugal_Force,
        Heat,
        Limit
    };
};
struct BoundaryTypeFlow {
    enum {
        Velo_X,
        Velo_Y,
        Velo_Z,
        Pressure,
        Temperature,
        Heat,
        Thermal_Flux,
        OtherScalar,
        Limit
    };
};
#endif	/* _BOUNDARYTYPE_H */
}
