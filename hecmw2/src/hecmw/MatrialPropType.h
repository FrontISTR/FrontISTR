/* 
 * File:   MatrialPropType.h
 * Author: ktakeda
 *
 * Created on 2009/05/14, 19:38
 */
#ifndef _MATRIAL_Prop_TYPE_H_317b7e89
#define	_MATRIAL_Prop_TYPE_H_317b7e89

namespace pmw{
struct MaterialPropType{
    enum{
        Density,     //ρ(密度)
        Poisson,     //ν(ポアソン比)
        YoungModule, //E(縦弾性係数)
        Temp_Depend_YoungModule,//温度依存 E(縦弾性係数)
        Temp_Depend_Poisson,    //温度依存 ν(ポアソン比)
        Linear_Expansion, //線膨張
        Thermal_Conductivity, //熱伝導
        Heat_Transfer_Rate,   //熱伝達
        Cp,//定圧比熱
        Cv //定積比熱
    };
};
}
#endif	/* _MATRIALTYPE_H */

