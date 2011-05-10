/* 
 * File:   BoundaryType.h
 * Author: ktakeda
 *
 * Modify     2009/05/25
 * Created on 2009/05/14, 14:13
 */

namespace pmw{
#ifndef _BOUNDARYTYPE_H_
#define	_BOUNDARYTYPE_H_
// 境界条件種類
struct BoundaryType{
    enum{
        Dirichlet,
        Neumann,
        Limit//終端
    };
};
// 境界メッシュのタイプ
struct BoundaryMeshType{
    enum{
        Node,
        Edge,
        Face,
        Volume,
        Limit//終端
    };
};

struct BoundaryTypeNode{
    enum {
        Load,
        Disp,
        Velo,
        Accel,
        Temp,
        Thermal_Flux,
        Limit//終端
    };
};
struct BoundaryTypeFace{
    enum {
        Pressure,
        TractionVector,
        Temp,
        Thermal_Flux,
        Limit//終端
    };
};
struct BoundaryTypeVolume{
    enum {
        Accel,
        Gravity,
        Centrifugal_Force,//遠心力
        Heat,
        Limit//終端
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
        OtherScalar,
        Limit//終端
    };
};
#endif	/* _BOUNDARYTYPE_H */
}


